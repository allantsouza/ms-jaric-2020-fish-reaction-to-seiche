# Detection of shifts
###CREATION OF THE MAIN DATA.TABLE - "shift_dataset"
# Tilt of thermocline calculation - order dataset on location_order in order to know to which direction the thermocline was tilted

thermocline <- merge(thermocline,
                     logger_positions[,.(location, lake, location_order)],
                     by = c("location", "lake"))

setkey(thermocline, lake, interval, location_order)

shift_dataset <- thermocline[,.(tilt = depth[1] - depth[2]),
                             by = .(lake, interval, step_order, slope, therm_part)]

#SETKEY, ensuring that the data table is sorted based on time
setkey(shift_dataset, lake, step_order, slope, therm_part, interval)

#Smoothing the tilt over time (2h moving window) #NOTE :
shift_dataset[, tilt_smooth := zoo::rollapply(width = 12,
                                              FUN = mean,
                                              data = tilt,
                                              partial = TRUE),
              by = .(lake, step_order, slope, therm_part)]

# Shift of thermocline calculation:
# shift (shift_smooth) is a distance of Seiche depth change between two consecutive records
# speed (shift_speed_smooth) is equal to height divided by the actual time passed between two consecutive records
shift_dataset[, shift_smooth := c(diff(tilt_smooth), NA),
              by = .(lake, step_order, slope, therm_part)]

shift_dataset[, shift_speed_smooth := (shift_smooth/as.numeric(difftime(interval[1],interval[2], units= c("hours")))),
              by = .(lake, step_order, slope, therm_part)]

#Identification of fast Seiche shift periods - "shift speed" has to be larger than a threshold speed (PAR_SHIFT_HEIGHT, which is X meters per hour)
shift_dataset[shift_speed_smooth <= PAR_SHIFT_SPEED & shift_speed_smooth >= -PAR_SHIFT_SPEED, shift_category := 0]
shift_dataset[shift_speed_smooth > PAR_SHIFT_SPEED, shift_category := 1]
shift_dataset[shift_speed_smooth < -PAR_SHIFT_SPEED, shift_category := -1]

#add elevating side - Since the dataset where ordered by location_order before computing the depth differences - assign order of logger line where it is elevating
shift_dataset[shift_speed_smooth > PAR_SHIFT_SPEED, pos_order_elevating := 2]
shift_dataset[shift_speed_smooth < PAR_SHIFT_SPEED, pos_order_elevating := 1]

#Numeration of each individual shift within groups
shift_dataset[, shift_phase := rleid(shift_category), by = .(lake, step_order, slope, therm_part)]

###CREATION OF SECOND DATA.TABLE WITH EXTRACTED DATA ON SHIFTS - "shifts_table"
#shift_dataset of waves (creation of "tilts and shifts" datatable)
shifts_table <- shift_dataset[,.(shift_duration = as.numeric(difftime(max(interval), min(interval), units=c("secs"))), 
                                 shift_height = abs(tilt_smooth[which.max(interval)] - tilt_smooth[which.min(interval)]), 
                                 max_shift_speed = max(abs(shift_speed_smooth)),
                                 max_shift_speed_time = interval[which.max(abs(shift_speed_smooth))], 
                                 start_time = min(interval),
                                 end_time = max(interval)
), by = .(lake, pos_order_elevating, step_order, slope, therm_part, shift_phase, shift_category)]

#Identification of the dataset of shifts (exclusion from the datatable of periods that are not recognized as shifts, keeping only the deepest thermocline for periods where more than one is identified)
shifts_table <- shifts_table[(shifts_table$shift_category == 1 | shifts_table$shift_category == -1) & step_order == 1]
#Exclusion of all shifts that were smaller than one meter of depth
shifts_table <- shifts_table[abs(shift_height)>=1]

#add elevating_side column
sides <- unique(thermocline[,.(lake, elevating_side = location,  pos_order_elevating = location_order)])
shifts.sides <- merge(sides, shifts_table, by = c("lake", "pos_order_elevating"), all.y = T)

ggplot(data = shifts_table,
       mapping = aes(x =shift_duration,
                     y = shift_height,
                     col = therm_part)) +
  geom_point() + 
  facet_wrap(~lake, scales = c("free"), ncol = 2)


#add unique shift id to table
shifts.sides[, shiftid := 1:.N]



ggplot(data = shift_dataset[therm_part == "center" &  step_order == 1 & lake == "Chabarovice" ])+
  geom_line(mapping = aes(x =interval , y = tilt))+
  geom_line(mapping = aes(x =interval , y = tilt_smooth, col = as.factor(shift_category), group = shift_phase))






