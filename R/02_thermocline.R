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
                     group = location)) +
  geom_line() +
  facet_wrap(~ month, scales = "free") +
  ylim(c(0, 10))

#Check errors - TODO:
aa <-  thermocline_temperatures_rolled_sub[!is.na(temperature) & slope == 1 & step_order == 1 & location %in% c("East", "West") & therm_part == "tstart"]
aa[month == 6 & depth < 2]
bb <- temperatures[location == "West" & ts == "2015-06-25 04:10:00"]
plot(bb$depth, bb$temperature)
temperatures[temperature == 18.2360]
ggplot(bb[location == "West" & ts == "2015-06-25 04:10:00"], aes(x = depth, y = temperature))+
  geom_point() +
  geom_hline(data = data.table(itp = c(18.2360)), aes(yintercept = itp)) +
  xlim(c(0, 15))

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
#get mean depth of thremocline in each 5min interval
therm_lake <- thermocline_location[,.(lake_therm_depth = mean(depth)), by = .(lake, interval,step_order, slope, therm_part)]
# smooth mean depth by therm_bal_ws days moving window
setkey(therm_lake,  interval)
therm_lake[, balanced_therm_depth := roll_time_window(span = PAR_THERMOCLINE_BALACE, FUN = mean, x = lake_therm_depth, times = interval),
           by = .(lake, step_order, slope, therm_part)]
therm_lake[, "lake_therm_depth" := NULL]

#therm mean now contains depth of thermocline in balanced state without wind
#compute thermocline thickness in each 5 minute interval
thermocline_location_wide <- dcast(data = thermocline_location, lake + location + interval + step_order ~ therm_part, value.var = "depth")
thermocline_location_wide[, thickness := end - start]

#merge data from locations to data from whole lake lake  - get deviation in locations of logger lines
th_mt <- merge(thermocline_location, therm_lake[therm_part == "center"],
               by = c("lake", "interval",  "step_order", "slope", "therm_part"))

th_mt[, deviation := balanced_therm_depth - depth]

th_m <- merge(th_mt, thermocline_location_wide, by = c("lake", "interval", "location", "step_order"))
th_m[, Month := month(interval)]

# Overview
ggplot(th_m[abs(deviation) > 1], aes(x = deviation, y = thickness, col = location)) +
  geom_point() +
  facet_wrap(~Month) + 
  geom_smooth(se = F, method = "lm")


ggplot(th_m[therm_part == "center"], aes(x = balanced_therm_depth, y = depth, col = location)) +
  geom_point() +
  facet_wrap(~location)

ggplot(th_m[therm_part == "center"], aes(x = interval, y = thickness, col = location)) + 
  geom_line() +
  facet_wrap(~Month, scales = "free", ncol = 1)

ggplot(th_m[therm_part == "center" ], aes(x = interval, y = deviation, col = location)) + 
  geom_line() +
  facet_wrap(~Month, scales = "free", ncol = 1)

therm_temp <- merge(thermocline_wide, th_m, by = c("location", "interval"))

ggplot(th_m[], aes(x = interval, y = balanced_therm_depth, col = therm_part))+
  geom_line()
