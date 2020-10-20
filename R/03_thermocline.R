# Script computing thermocline



# Finding thermocline -----------------------------------------------------
#hobo_data[, rho := 999,974950 *(1 - (temperature+288.9414)/(508929.2*(temperature+68.12963))*(temperature-3.9863)^2)]



hobo_data <- data.table(load_hobo_data())
temperatures_monotonic <- data.table(load_temperature_data())

a1 <- -3.983035
a2 <- 301.797
a3 <- 522528.9
a4 <- 69.34881
a5 <- 999.974950
hobo_data[, water_density := a5 * (1 - ((((a1 + temperature)^2) * (a2 + temperature))/(a3*(temperature+a4))))] #http://metgen.pagesperso-orange.fr/metrologieen19.htm

#remove intervals were there is no stratification (difference in temperature on surface and bottom < 2) or there is not enough logs (<13)
setkey(hobo_data, location, interval, depth)
hobo_data[, out := .N < 13 | mean(temperature[2:4])- mean(tail(temperature, 2)) < 2, by =.(location, interval)]
hobo_data_therm <- hobo_data[out == F]



# Compute thermocline -----------------------------------------------------

thermocline <- hobo_data_therm[, compute_thermocline(depth = depth,
                                                     temperature = temperature,
                                                     diff_threshold = PAR_THERMOCLINE_SLOPE,
                                                     depth_res = 0.1),
                               by = .(thermocline_ts = interval, lake, location)]
thermocline$slope <- PAR_THERMOCLINE_SLOPE

ggplot(data = thermocline[step_order == 1 & location %in% c("East", "West") & thermocline_ts %between% c("2015-07-13", "2015-07-18")],
       mapping = aes(x = thermocline_ts,
                     y = temperature_start,
                     col = location)) +
  geom_point(shape = ".") +
  ylab("Temperature") +
  xlab("Date") +
  ggtitle("Top of the thermocline over time")


# Investigation why some thermocline starts at 12 Celsius
wierd_interval <- thermocline[temperature_start < 12 & step_order == 1 & location %in% c("East") & thermocline_ts %between% c("2015-07-13", "2015-07-18"),][6]
wierd_profile <- hobo_data_therm[interval == wierd_interval$thermocline_ts & location == wierd_interval$location]
wierd_profile_smooth <- data.frame(interpolate_temperature_profile(depth = wierd_profile$depth, temperature = wierd_profile$temperature))
wierd_profile_smooth$is_thermocline <- wierd_profile_smooth$slope < -PAR_THERMOCLINE_SLOPE
ggplot(data = wierd_profile,
       mapping = aes(x = temperature,
                     y = depth)) +
  scale_y_reverse()+
  geom_point() +
  geom_line(data =  wierd_profile_smooth) + 
  geom_point(data = wierd_profile_smooth[], aes(group = rleid(is_thermocline), col = is_thermocline), alpha=0.3)



# Seasonal thermocline ----------------------------------------------------


# calculation of mean thermocline temperature
setkey(thermocline, location, step_order, slope, thermocline_ts)
window_size <- 3*86400
# tcenter
# thermocline[, tcenter := roll_time_window(x = (temperature_start+temperature_end)/2, 
#                                        times = thermocline_ts,
#                                        span = window_size,
#                                        FUN = function(x) median(x, na.rm = T)),
#             by = .(location, step_order, slope)]
# 
# 

# tstart
thermocline[, tstart_smoothed := roll_time_window(x = temperature_start,
                                       times = thermocline_ts,
                                       span = window_size,
                                       FUN = function(x) median(x, na.rm = T)),
                        by = .(location, step_order, slope)]

# tend
thermocline[, tend_smoothed := roll_time_window(x = temperature_end, 
                                     times = thermocline_ts,
                                     span = window_size,
                                     FUN = function(x) median(x, na.rm = T)),
            by = .(location, step_order, slope)]
# tcrit
thermocline[, tcrit_smoothed := roll_time_window(x = temperature_crit, 
                                      times = thermocline_ts,
                                      span = window_size,
                                      FUN = function(x) median(x, na.rm = T)),
            by = .(location, step_order, slope)]

thermocline[, tcenter_smoothed := (tend_smoothed + tstart_smoothed)/2]

write_csv(thermocline, file = here("data","products", "thermocline_slope.csv"))

# melt so the roll is in long format
ggplot(data = thermocline[step_order == 1 & location %in% c("East", "West") ],
       mapping = aes(x = thermocline_ts,
                     y = tcenter_smoothed,
                     col = location)) +
  geom_point(shape = ".")

roll_cols <- c("tcenter_smoothed", "tstart_smoothed", "tend_smoothed", "tcrit_smoothed")
thermocline_full <- melt(thermocline,
                         id.vars = c("lake", "location", "thermocline_ts", "step_order", "slope"),
                         variable.name = "therm_part",
                         value.name = "temperature_smoothed",
                         measure.vars = roll_cols)

thermocline_full[, therm_part := gsub("_smoothed", "", x = therm_part)]

#TODO: overview of gaps in seconds - should be 0! or interpolate otherwise to have full dataset
thermocline_full[step_order == 1,
                 .(thermocline_ts, time_difference = c(0, diff(as.numeric(thermocline_ts)))),
                 by = .(lake, location, step_order, slope, therm_part)][time_difference > 300]


# TODO: this is very important step! Getting one temperature for the whole lake!
# calculating means of the two locations for each time
thermocline_full[, lake_therm_temperature_smoothed := mean(temperature_smoothed), by = .(lake, thermocline_ts, step_order, slope, therm_part)]
#exclusion of periods with no thermocline
thermocline_full <- thermocline_full[!is.na(temperature_smoothed)]
ggplot(data = thermocline_full[step_order == 1 & location %in% c("East", "West")],
       mapping = aes(x = thermocline_ts,
                     y = temperature_smoothed,
                     col = therm_part,
                     linetype  = location)) +
  geom_line() +
  xlab("Date") +
  ylab("Temperature (moving median)")



# Get depth of thermocline ------------------------------------------------
setDT(temperatures_monotonic)
temperatures_monotonic_strictly <- temperatures_monotonic[, .(depth = max(depth)), by = .(ts, location, temperature)]
setkey(temperatures_monotonic_strictly, location, ts, depth)
temperatures_monotonic_strictly$slope <- NULL
#temperatures_monotonic_strictly <- temperatures_monotonic_strictly[depth > 3]
temperatures_monotonic_strictly[, temperature_strictly_decreasing := temperature - 0.001 * (1:.N), by = .(ts, location)]
# Previous code computed lake-wide temperature of thermocline for each 5min interval
# Now get depth of occurence of that temperature
thermocline_full[, ':=' (temperature_roll = lake_therm_temperature_smoothed)]
temperatures_monotonic_strictly[, ':=' (thermocline_ts = ts, temperature_roll = temperature_strictly_decreasing)]
setkey(thermocline_full, location, thermocline_ts, temperature_roll)
setkey(temperatures_monotonic_strictly, location, thermocline_ts, temperature_roll)
thermocline_temperatures_rolled <- temperatures_monotonic_strictly[thermocline_full, , on = c("location", "thermocline_ts", "temperature_roll"), roll = "nearest"]
thermocline_temperatures_rolled[, temperature_roll := NULL]

# Hotfix cure for some wierd profiles 
thermocline_temperatures_rolled[depth > 15, depth := 15]

#remove joins futher that 15 minutes
thermocline_temperatures_rolled <- thermocline_temperatures_rolled[abs(as.numeric(ts) - as.numeric(thermocline_ts)) < 60*15 ]

#Thermocline dynamics - PLOT 
ggplot(data = thermocline_temperatures_rolled[!is.na(lake_therm_temperature_smoothed) & step_order == 1 ],
       mapping = aes(x = thermocline_ts,
                     y = depth,
                     col = therm_part)) +
  geom_line() +
  facet_wrap(~ location, ncol = 1) +
  scale_y_reverse(limits=c(20,0)) +
  xlab("Date") + 
  ylab("Depth") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, name = "Spectral"),
                     breaks = c("tstart", "tcrit", "tcenter", "tend"),
                     labels = c("start", "crit", "center", "end")) +
  labs(col = "Thermocline part")

#extract columns of interest
thermocline_location <- thermocline_temperatures_rolled[step_order == 1, .(lake,
                                                             location,
                                                             thermocline_ts,
                                                             depth,
                                                             temperature,
                                                             step_order,
                                                             lake_therm_temperature_smoothed,
                                                             slope,
                                                             therm_part = sub(therm_part, pattern = "^t(.*)$", replacement = "\\1"))]


# Calculate mean seasonal thermocline for whole lake
# calculation of mean thermocline temperature for every single day (mean from both lines, calculated as diffrence between up_depth and down_depth)
# Get mean depth of thermocline in the whole lake
setkey(thermocline_location, location, thermocline_ts)
# Get mean depth of thremocline in each 5min interval
therm_lake <- thermocline_location[,.(lake_therm_depth = mean(depth)), by = .(lake, thermocline_ts, step_order, slope, therm_part, lake_therm_temperature_smoothed)]
# smooth mean depth by PAR_THERMOCLINE_SMOOTH days moving window
setkey(therm_lake,  thermocline_ts)
therm_lake[, lake_therm_depth_smoothed := roll_time_window(span = PAR_THERMOCLINE_SMOOTH, FUN = mean, x = lake_therm_depth, times = thermocline_ts),
           by = .(lake, step_order, slope, therm_part)]
therm_lake[, "lake_therm_depth" := NULL]

ggplot(therm_lake[step_order == 1], aes(x = thermocline_ts, y = lake_therm_depth_smoothed, col = therm_part))+
  geom_point(shape = ".") +
  xlab("Date") +
  ylab("Depth") + 
  scale_y_reverse() +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, name = "Spectral"),
                     breaks = c("start", "crit", "center", "end"),
                     labels = c("start", "crit", "center", "end"))


# Get deviation from balanced depth
th_deviation <- merge(thermocline_location, therm_lake[, .(lake, thermocline_ts, slope, step_order, therm_part, lake_therm_depth_smoothed)],
               by = c("lake", "thermocline_ts",  "step_order", "slope", "therm_part"))
th_deviation[, deviation := lake_therm_depth_smoothed - depth]


# Compute thermocline thickness
thermocline_location_wide_depth <- dcast(data = thermocline_location, lake + location + slope + thermocline_ts + step_order ~ therm_part, value.var = "depth")
thermocline_location_wide_depth[, thickness := end - start]

thermocline_data_thickness <- merge(th_deviation, thermocline_location_wide_depth[, .(lake, thermocline_ts, location, step_order, slope, thickness)],
                          by = c("lake", "thermocline_ts", "location", "step_order", "slope"))

# Compute thermocline strength
thermocline_location_wide_temperature <- dcast(data = thermocline_location, lake + location + slope + thermocline_ts + step_order ~ therm_part, value.var = "temperature")
thermocline_location_wide_temperature[, strength := start - end]

thermocline_data <- merge(thermocline_data_thickness, thermocline_location_wide_temperature[, .(lake, thermocline_ts, location, step_order, slope, strength)],
                          by = c("lake", "thermocline_ts", "location", "step_order", "slope"))



write_csv(x = thermocline_data %>% filter(thermocline_ts %between% DATE_RANGE), file = here("data", "products", "thermocline_data.csv"))

# Overview

# Deviation x Thickness
ggplot(thermocline_data[abs(deviation) > 0.7], aes(x = deviation, y = thickness, col = location)) +
  geom_point(pch = ".") + 
  geom_smooth(method = "lm")

# Lake therm depth vs instant depth (basically deviation)
ggplot(thermocline_data[therm_part == "center"], aes(x = lake_therm_depth_smoothed, y = depth, col = location)) +
  geom_point(pch = ".") + 
  facet_wrap(~location)

# Overview of thickness
ggplot(thermocline_data[therm_part == "center",], aes(x = thermocline_ts, y = thickness, col = location)) + 
  geom_line() +
  geom_vline(data = unique(thermocline_data[therm_part == "center", .(as.POSIXct(as.Date(thermocline_ts)))]), aes(xintercept = V1))

# Overview of deviation
ggplot(thermocline_data[therm_part == "center" ], aes(x = thermocline_ts, y = deviation, col = location)) + 
geom_line()

# TODO: deviation and thickness has different frequency!

# Overview of strength
ggplot(thermocline_data[therm_part == "center",], aes(x = thermocline_ts, y = strength, col = location)) + 
  geom_line()

# Strengh vs deviation - not correlated - good ;)
ggplot(thermocline_data[therm_part == "center",], aes(x = deviation, y = strength, col = location)) + 
  geom_point(pch = ".") +
  geom_density2d()

# Thickness vs deviation - not correlated - good ;)
ggplot(thermocline_data[therm_part == "center",], aes(x = deviation, y = thickness, col = location)) + 
  geom_point(pch = ".") +
  geom_density2d()
