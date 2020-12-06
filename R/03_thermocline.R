# Script computing thermocline



# Computation of thermocline for each time interval -----------------------------------------------------

temperatures_monotonic <- data.table(load_temperature_data())

thermocline <- temperatures_monotonic %>%
  #filter(ts %between% c("2015-07-02", "2015-07-29")) %>%
  group_nest(location, ts) %>%
  mutate(thermocline = map(.x = data, .f = safely(compute_thermocline_thermo_depth, otherwise = NULL))) %>%
  mutate(thermocline_layer_result = map(.x = thermocline, ~ .x$result)) %>%
  mutate(stratification = map(.x = data, .f = function(x) sd(x$temperature[which(x$depth < 17.5)]))) %>% # measure of stratification
  dplyr::select(-data, -thermocline) %>%
  unnest(c("thermocline_layer_result", "stratification")) %>%
  write_csv(here("data", "products", "thermocline_thermo_depth.csv"))

#thermocline <- read_csv(here("data", "products", "thermocline_wtr_layer.csv"))
#thermocline <- read_csv(here("data/products/thermocline_themo_depth.csv"))

ggplot(data = thermocline %>% filter(ts %between% c("2015-08-15", "2015-08-27")),
       mapping = aes(x = ts,
                     y = temperature_start,
                     col = location)) +
  geom_point(shape = ".") +
  geom_density2d()+
  ylab("Temperature") +
  xlab("Date") +
  ggtitle("Top of the thermocline over time")

# Smoothing of thermocline in time ----------------------------------------------------
# Since the method is not very stable on such a various number of different profiles
# the thermocline is smoothed by finding isotherms in moving window

# Calculation of median thermocline temperature (moving window)
window_size <- 3*86400

thermocline_smoothed <- thermocline %>% 
  pivot_longer(cols = matches("depth|temperature"),
                   names_to = c("var", "therm_part"), 
                   names_pattern = "(.*)_(.*)")  %>%
  pivot_wider(id_cols = c("location", "ts", "therm_part"), names_from = "var") %>%
  group_by(location, therm_part) %>%
  arrange(ts) %>%
  mutate(temperature_smoothed = roll_time_window(x = temperature,
                          times = ts,
                          span = window_size,
                          FUN = function(x) median(x, na.rm = T))) %>%
  ungroup()

# Preview
ggplot(data = thermocline_smoothed,
       mapping = aes(x = ts,
                     y = temperature,
                     col = therm_part)) +
  geom_point(shape = ".") +
  geom_line(aes(y = temperature_smoothed)) +
  facet_wrap(~ location, ncol = 1)



# Get depth of thermocline based on smoothed temperature 
setDT(temperatures_monotonic)
temperatures_monotonic_strictly <- temperatures_monotonic[, .(depth = max(depth)), by = .(ts, location, temperature)]
setkey(temperatures_monotonic_strictly, location, ts, depth)

# Make profile strictly monotonic
temperatures_monotonic_strictly[, temperature_strictly_decreasing := temperature - 0.001 * (1:.N), by = .(ts, location)]

# This is very important step!
# Getting one temperature for the whole lake! (computed as mean)
thermocline_full <- as.data.table(thermocline_smoothed)
#thermocline_full[, lake_therm_temperature_smoothed := mean(temperature_smoothed), by = .(ts, therm_part)]
thermocline_full[, lake_therm_temperature_smoothed := temperature_smoothed, by = .(ts, therm_part)]

thermocline_full[, ':=' (temperature_roll = lake_therm_temperature_smoothed,
                         thermocline_ts = ts,
                         depth = NULL, 
                         temperature = NULL, 
                         ts = NULL)]

temperatures_monotonic_strictly[, ':=' (thermocline_ts = ts, temperature_roll = temperature_strictly_decreasing)]
setkey(thermocline_full, location, thermocline_ts, temperature_roll)
setkey(temperatures_monotonic_strictly, location, thermocline_ts, temperature_roll)
thermocline_temperatures_rolled <- temperatures_monotonic_strictly[thermocline_full, , on = c("location", "thermocline_ts", "temperature_roll"), roll = "nearest"]
thermocline_temperatures_rolled[, temperature_roll := NULL]

# Hotfix cure for some wierd profiles 
thermocline_temperatures_rolled[depth > 15, depth := 15]
thermocline_temperatures_rolled <- thermocline_temperatures_rolled[depth > 4]


# Remove joins futher that 15 minutes
#thermocline_temperatures_rolled <- thermocline_temperatures_rolled[abs(as.numeric(ts) - as.numeric(thermocline_ts)) < 60*15 ]

# Thermocline dynamics - PLOT 
ggplot(data = thermocline_temperatures_rolled[],
       mapping = aes(x = thermocline_ts,
                     y = depth,
                     col = therm_part)) +
  geom_line() +
  facet_wrap(~ location, ncol = 1) +
  scale_y_reverse(limits=c(20,0)) +
  xlab("Date") + 
  ylab("Depth") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, name = "Spectral"),
                     breaks = c("start", "crit", "end"),
                     labels = c("start", "crit", "end")) +
  labs(col = "Thermocline part")

# Extract columns of interest
thermocline_location <- thermocline_temperatures_rolled[, .(location,
                                                             thermocline_ts,
                                                             depth,
                                                             temperature,
                                                             therm_part)]


# Compute derived thermocline attributes ----------------------------------
# Approximate linearly missing records - not a big deal - few hours
therm_lake <- thermocline_location %>%
  as_tibble() %>%
  mutate(ts_num = as.numeric(thermocline_ts) - min(as.numeric(thermocline_ts))) %>%
  group_nest(therm_part, location) %>%
  mutate(full_time_series = map(.x = data, .f = function(x){
   tibble(
     ts_num = seq(0, max(x$ts_num), 300),
     depth = approx(x = x$ts_num, y = x$depth, xout = seq(0, max(x$ts_num), 300))$y,
     temperature = approx(x = x$ts_num, y = x$temperature, xout = seq(0, max(x$ts_num), 300))$y,
     thermocline_ts = min(x$thermocline_ts) + seq(0, max(x$ts_num), 300)
   )
   
  })) %>%
  select(-data) %>%
  unnest(c("full_time_series"))

# Calculate seasonal thermocline and deviation
therm_lake_seasonal <-  therm_lake %>%
  as_tibble() %>%
  group_nest(therm_part, location) %>%
  mutate(location_therm_depth_smoothed = map(.x = data, .f = function(x){
    gam_model <- mgcv::gam(depth ~ s(ts_num, k = 100), data = x)
    predict(gam_model, newdata = data.frame(ts_num = x$ts_num))
  })) %>%
  unnest(c("data", "location_therm_depth_smoothed"))

# Smooth deviation by fft
# Function to smooth variable using FFT
# returns list of fft result (freq - frequency in hours, spec - frequency spectrum) and smoothed variable
smooth_fft <- function(x, cutoff = 700) {
  x_fft <- fft(x)
  x_cutoffed <- x_fft
  x_cutoffed[(cutoff):(length(x)-cutoff)] <- 0 + 0i 
  tibble(fft_freq = 1/((1:length(x)) * 12/length(x)), 
         fft_spec = Re(x_fft),
         fft_smooth = Re(fft(x_cutoffed, inverse = TRUE)/length(x_fft)))
}


therm_lake_deviation_t <- therm_lake_seasonal %>%
  mutate(deviation = depth - location_therm_depth_smoothed) %>%
  arrange(location, therm_part, thermocline_ts) %>%
  group_nest(location, therm_part) %>%
  mutate(deviation_smoothed = map(.x = data, .f = ~ smooth_fft(.x$deviation))) %>%
  unnest(c("data", "deviation_smoothed")) %>%
  rename(deviation_fft_smoothed = fft_smooth) %>%
  mutate(depth_fft_smoothed = location_therm_depth_smoothed + deviation_fft_smoothed)

therm_lake_deviation_t %>% 
  select(location, therm_part, fft_freq, fft_spec) %>%
  write_csv("data/products/thermocline_fft.csv")

ggplot(therm_lake_deviation_t, 
       aes(x = thermocline_ts, y = depth, col = therm_part)) +
  geom_line(alpha = 0.5) +
  geom_line(aes(y = location_therm_depth_smoothed), linetype = 2) +
  geom_line(aes(y = depth_fft_smoothed)) +
  xlab("Date") +
  ylab("Depth") + 
  scale_y_reverse() +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, name = "Spectral"),
                     breaks = c("start", "crit" , "end"),
                     labels = c("start", "crit", "end")) +
  facet_wrap(~ location, ncol = 1)


# Prepare temperatures for merging

therm_lake_deviation_tt <- therm_lake_deviation_t %>% mutate(depth_fft_smoothed_rounded = round(depth_fft_smoothed, 1))

therm_lake_deviation <- temperatures_monotonic_strictly %>% as_tibble() %>%
  mutate(depth_fft_smoothed_rounded = round(depth, 1)) %>%
  select(ts, location, depth, depth_fft_smoothed_rounded, temperature_strictly_decreasing) %>%
  rename(thermocline_ts = ts,
         temperature_fft_smoothed = temperature_strictly_decreasing) %>%
  right_join(therm_lake_deviation_tt, by = c("location", "thermocline_ts", "depth_fft_smoothed_rounded")) %>%
  select(-depth_fft_smoothed_rounded)


therm_lake_deviation %>% filter(thermocline_ts == "2015-06-09 00:15:00" & depth_fft_smoothed_rounded == 6.8) 
b$depth_fft_smoothed_rounded 

ggplot(therm_lake_deviation, 
       aes(x = thermocline_ts, y = temperature, col = therm_part)) +
  geom_line(alpha = 0.5) +
  geom_line(aes(y = temperature_fft_smoothed)) +
  xlab("Date") +
  ylab("Temperature") +
  facet_wrap(~ location, ncol = 1) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, name = "Spectral"),
                     breaks = c("start", "crit" , "end"),
                     labels = c("start", "crit", "end")) 


# thickness - difference of depth between start and end
# stregth - difference of degrees between start and end
# mean_gradient - temperature decrease in degrees/m

thermocline_lake_thickness_strength <- therm_lake_deviation %>% 
  pivot_wider(id_cols = c("thermocline_ts", "location"),
              names_from = "therm_part",
              values_from = c("temperature", "depth", "depth_fft_smoothed", "temperature_fft_smoothed")) %>%
  mutate(thickness = depth_end - depth_start,
         thickness_fft_smoothed = depth_fft_smoothed_end - depth_fft_smoothed_start,
         strength = temperature_start - temperature_end,
         strength_fft_smoothed = temperature_fft_smoothed_start - temperature_fft_smoothed_end) %>%
  select(thermocline_ts, location, thickness, thickness_fft_smoothed, strength) %>%
  inner_join(therm_lake_deviation, by = c("location", "thermocline_ts")) %>%
  mutate(mean_gradient = strength/thickness, 
         mean_gradient_fft_smoothed = strength_fft_smoothed / thickness_fft_smoothed)


# Export
write_csv(x = thermocline_lake_thickness_strength %>% filter(thermocline_ts %between% DATE_RANGE & !is.na(mean_gradient)),
          file = here("data", "products", "thermocline_data.csv"))
