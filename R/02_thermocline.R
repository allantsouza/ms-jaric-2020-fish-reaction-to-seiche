# Script computing thermocline

#finding thermocline ####
#hobo_data[, rho := 999,974950 *(1 - (temperature+288.9414)/(508929.2*(temperature+68.12963))*(temperature-3.9863)^2)]
a1 <- -3.983035
a2 <- 301.797
a3 <- 522528.9
a4 <- 69.34881
a5 <- 999.974950
hobo_data[, rho := a5 * (1 - ((((a1 + temperature)^2) * (a2 + temperature))/(a3*(temperature+a4))))] #http://metgen.pagesperso-orange.fr/metrologieen19.htm

#remove intervals were there is no stratification (difference in temperature on surface and bottom < 2) or there is not enough logs (<13)
setkey(hobo_data, location, interval, depth)
hobo_data[, out := .N < 13 | mean(temperature[2:4])- mean(tail(temperature, 2)) < 2, by =.(location, interval)]
hobo_data_therm <- hobo_data[out == F]

# compute thermocline -------------------------------------------------------------------------------------
thermocline <- hobo_data_therm[, compute_thermocline(depth = depth,
                                                     temperature = temperature,
                                                     diff_threshold = PAR_THERMOCLINE_SLOPE,
                                                     depth_res = 0.1),
                               by = .(interval, lake, location)]
thermocline$slope <- PAR_THERMOCLINE_SLOPE

# TODO: export this dataset as thermocline_5min

# Get seasonal thermocline --------------------------------------------------------------------------------

# calculation of mean thermocline temperature for every single day (mean from both lines, calculated as diffrence between up_depth and down_depth)
setkey(thermocline, location, step_order, slope, interval)
window_size <- 14*86400
# tcenter
thermocline[, tcenter := roll_time_window(x = (temperature_start+temperature_end)/2, 
                                       times = interval,
                                       span = window_size,
                                       FUN = function(x) median(x, na.rm = T)),
            by = .(location, step_order, slope)]

# tstart
thermocline[, tstart := roll_time_window(x = temperature_start,
                                       times = interval,
                                       span = window_size,
                                       FUN = function(x) median(x, na.rm = T)),
                        by = .(location, step_order, slope)]

# tend
thermocline[, tend := roll_time_window(x = temperature_end, 
                                     times = interval,
                                     span = window_size,
                                     FUN = function(x) median(x, na.rm = T)),
            by = .(location, step_order, slope)]
# tcrit
thermocline[, tcrit := roll_time_window(x = temperature_crit, 
                                      times = interval,
                                      span = window_size,
                                      FUN = function(x) median(x, na.rm = T)),
            by = .(location, step_order, slope)]


#melt so the roll is in long format
ggplot(data = thermocline[step_order == 1 & location %in% c("East", "West") ],
       mapping = aes(x = interval,
                     y = tstart,
                     col = location)) +
  geom_point(shape = ".")

roll_cols <- c("tcenter", "tstart", "tend", "tcrit")
thermocline_full <- melt(thermocline,
                         id.vars = c("lake", "location", "interval", "step_order", "slope"),
                         variable.name = "therm_part",
                         value.name = "temperature_smooth",
                         measure.vars = roll_cols)

#TODO: overview of gaps in seconds - should be 0! or interpolate otherwise to have full dataset
thermocline_full[step_order == 1,
                 .(interval, time_difference = c(0, diff(as.numeric(interval)))),
                 by = .(lake, location, step_order, slope, therm_part)][time_difference > 300]


# TODO: this is very important step! Getting one temperature for the whole lake!
# calculating means of the two locations for each time
thermocline_full[, temperature_lake_mean := mean(temperature_smooth), by = .(lake, interval, step_order, slope, therm_part)]
#exclusion of periods with no thermocline
thermocline_full <- thermocline_full[!is.na(temperature_smooth)]
ggplot(data = thermocline_full[step_order == 1 & location %in% c("East", "West")],
       mapping = aes(x = interval,
                     y = temperature_smooth,
                     col = therm_part,
                     linetype  = location))+
  facet_wrap(~ slope, ncol = 1)


## Get depth for thermocline ####
# Previous code computed lake-wide temperature of thermocline for each 5min interval
# Now get depth of occurence of that temperature
thermocline_full[, ':=' (temperature = temperature_lake_mean)]
temperatures_monotonic[, ':=' (interval = ts)]
setkey(thermocline_full, location, interval, temperature)
setkey(temperatures_monotonic, location, interval, temperature)
thermocline_temperatures_rolled <- temperatures_monotonic[thermocline_full,, on = c("location", "interval", "temperature"), roll = "nearest"]
#remove joins futher that 15 minutes
thermocline_temperatures_rolled <- thermocline_temperatures_rolled[abs(as.numeric(ts) - as.numeric(interval)) < 60*15 ]

#Thermocline dynamics - PLOT 
ggplot(data = thermocline_temperatures_rolled[!is.na(temperature) & step_order == 1 ],
       mapping = aes(x = interval,
                     y = depth,
                     col = therm_part)) +
  geom_line() +
  facet_wrap(~ location, ncol = 1) +
  theme_minimal() + 
  scale_y_reverse() +
  xlab("Date") + 
  ylab("Depth") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, name = "Spectral"),
                     breaks = c("tstart", "tcrit", "tcenter", "tend"),
                     labels = c("start", "crit", "center", "end")) +
  labs(col = "Thermocline part")

#extract columns of interest
thermocline_location <- thermocline_temperatures_rolled[step_order == 1, .(lake,
                                                             location,
                                                             interval,
                                                             depth,
                                                             step_order,
                                                             temperature,
                                                             slope,
                                                             therm_part = sub(therm_part, pattern = "^t(.*)$", replacement = "\\1"))]


# Calculate mean seasonal thermocline for whole lake
# calculation of mean thermocline temperature for every single day (mean from both lines, calculated as diffrence between up_depth and down_depth)
# Get mean depth of thermocline in the whole lake
setkey(thermocline_location, location, interval)
# Get mean depth of thremocline in each 5min interval
therm_lake <- thermocline_location[,.(lake_therm_depth = mean(depth)), by = .(lake, interval, step_order, slope, therm_part)]
# smooth mean depth by therm_bal_ws days moving window
setkey(therm_lake,  interval)
therm_lake[, balanced_therm_depth := roll_time_window(span = PAR_THERMOCLINE_BALACE, FUN = mean, x = lake_therm_depth, times = interval),
           by = .(lake, step_order, slope, therm_part)]
therm_lake[, "lake_therm_depth" := NULL]

ggplot(therm_lake[step_order == 1], aes(x = interval, y = lake_therm_depth , col = therm_part))+
  geom_point(shape = ".") +
  geom_line(aes(y = balanced_therm_depth)) + 
  theme_minimal() + 
  xlab("Date") +
  ylab("Depth") + 
  scale_y_reverse() +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, name = "Spectral"),
                     breaks = c("start", "crit", "center", "end"),
                     labels = c("start", "crit", "center", "end"))


# Get deviation from balanced depth
th_deviation <- merge(thermocline_location, therm_lake[therm_part == "center", .(lake, interval, slope, step_order, balanced_therm_depth)],
               by = c("lake", "interval",  "step_order", "slope"))
th_deviation[, deviation := balanced_therm_depth - depth]


# Compute thermocline thickness
thermocline_location_wide <- dcast(data = thermocline_location, lake + location + slope + interval + step_order ~ therm_part, value.var = "depth")
thermocline_location_wide[, thickness := end - start]

thermocline_data <- merge(th_deviation, thermocline_location_wide[, .(lake, interval, location, step_order, slope, thickness)],
                          by = c("lake", "interval", "location", "step_order", "slope"))

write_csv(x = thermocline_data, path = here("data", "products", "thermocline_data.csv"))

# Overview
ggplot(thermocline_data[abs(deviation) > 0.7], aes(x = deviation, y = thickness, col = location)) +
  geom_point() +
  geom_smooth(se = F, method = "lm") +
  theme_minimal()

ggplot(thermocline_data[therm_part == "center"], aes(x = balanced_therm_depth, y = depth, col = location)) +
  geom_point() +
  facet_wrap(~location) +
  theme_minimal()

ggplot(thermocline_data[therm_part == "center"], aes(x = interval, y = thickness, col = location)) + 
  geom_line() +
  theme_minimal()

ggplot(thermocline_data[therm_part == "center" ], aes(x = interval, y = deviation, col = location)) + 
  geom_line() +
  theme_minimal()