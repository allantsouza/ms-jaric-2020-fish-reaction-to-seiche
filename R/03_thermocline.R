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

# Remove joins futher that 15 minutes
thermocline_temperatures_rolled <- thermocline_temperatures_rolled[abs(as.numeric(ts) - as.numeric(thermocline_ts)) < 60*15 ]

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

# Calculate seasonal thermocline
therm_lake <- thermocline_location %>%
  as_tibble() %>%
  mutate(ts_num = as.numeric(thermocline_ts) - min(as.numeric(thermocline_ts))) %>%
  group_nest(therm_part, location) %>%
  mutate(lake_therm_depth_smoothed = map(.x = data, .f = function(x){
    gam_model <- mgcv::gam(depth ~ s(ts_num, k = 100), data = x)
    predict(gam_model, newdata = data.frame(ts_num = x$ts_num))
  })) %>%
  unnest(c("data", "lake_therm_depth_smoothed")) %>%
  as.data.table()

print(therm_lake[, .N, by = .(location, therm_part)])

ggplot(therm_lake, aes(x = thermocline_ts, y = depth, col = therm_part))+
  geom_point(shape = ".") +
  geom_line(aes(y = lake_therm_depth_smoothed)) +
  xlab("Date") +
  ylab("Depth") + 
  scale_y_reverse() +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, name = "Spectral"),
                     breaks = c("start", "crit" , "end"),
                     labels = c("start", "crit", "end")) + facet_wrap(~ location, ncol = 1)


# deviation - difference of thermocline depth from seasonal thermocline depth
th_deviation <- merge(thermocline_location, unique(therm_lake[, .(thermocline_ts, location, therm_part, lake_therm_depth_smoothed)]),
               by = c("thermocline_ts","location", "therm_part"))
th_deviation[, deviation := depth - lake_therm_depth_smoothed]


# thickness - difference of depth between start and end
thermocline_location_wide_depth <- dcast(data = thermocline_location,
                                        formula = location +  thermocline_ts ~ therm_part,
                                        value.var = "depth")

thermocline_location_wide_depth[, thickness := end - start]

thermocline_data_thickness <- merge(th_deviation,
                                    thermocline_location_wide_depth[, .(thermocline_ts,
                                                                        location,
                                                                        thickness)],
                          by = c("thermocline_ts", "location"))

# stregth - difference of degrees between start and end
thermocline_location_wide_temperature <- dcast(data = thermocline_location, location + thermocline_ts ~ therm_part, value.var = "temperature")
thermocline_location_wide_temperature[, strength := start - end]

thermocline_data <- merge(thermocline_data_thickness, thermocline_location_wide_temperature[, .(thermocline_ts, location, strength)],
                          by = c("thermocline_ts", "location"))


# mean_gradient - temperature decrease in degrees/m
thermocline_data[, mean_gradient := strength/thickness]


# Export
write_csv(x = thermocline_data %>% filter(thermocline_ts %between% DATE_RANGE),
          file = here("data", "products", "thermocline_data.csv"))
