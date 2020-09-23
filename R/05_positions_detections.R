# Script for processing raw detections and positions and merging with thermocline
# Assigning rough lateral position (distance from Western shore) to each detection and merging with thermocline. 
# OUTPUT: detections with distance from shore and interpolated values of thermocline in their specific location


# Loading data ------------------------------------------------------------


fish_raw <- read_csv(file = "data/raw/fishIDs.csv", col_types = "ccd")
tag_sns <- fish_raw$tag_sn


detection_dir <- "data/raw/db/fish/detections/"
position_dir <- "data/raw/db/fish/positions/"


lake_axis <- here("data/raw/sp/lake_axis.shp") %>%
  st_read()

lake_width <- as.numeric(st_length(lake_axis))

point_zero <- here("data/raw/sp/point_zero.shp") %>%
  st_read()


logger_pos <- here("data/raw/db/logger_positions.shp") %>% st_read()
logger_pos$dist_from_zero <- compute_distance_from_point_zero(positions = logger_pos, lake_axis = lake_axis, point_zero = point_zero)
distance_west <- logger_pos[logger_pos$locatin == "West",]$dist_from_zero
distance_east <- logger_pos[logger_pos$locatin == "East",]$dist_from_zero


thermocline <- here("data/products/thermocline_data.csv") %>%
  read_csv() %>%
  filter(step_order == 1 & slope == PAR_THERMOCLINE_SLOPE)

thermocline_wide <- thermocline %>%
  pivot_wider(id_cols = c("interval", 
                          "therm_part",
                          "lake_therm_depth_smoothed",
                          "lake_therm_temperature_smoothed"), names_from = "location",
              values_from = c("deviation", "thickness", "temperature", "depth"))




# For each fish -----------------------------------------------------------


for(i in 1:length(tag_sns)){
  positions <- read_csv(file = here(position_dir, paste0(tag_sns[i], ".csv")), col_types = c(up_depth = "d")) %>% 
    st_as_sf(coords = c("x", "y"), crs = 32633)
  
  detections <- read_csv(file = here(detection_dir, paste0(tag_sns[i], ".csv")))
  
  if(nrow(positions) > 0){
    # get distance from shore
    positions$dist_from_zero <- compute_distance_from_point_zero(positions = positions,
                                                                 lake_axis = lake_axis,
                                                                 point_zero = point_zero)
    
    # rolljoin positions and detections
    positions_dt <- as.data.table(positions)
    detections_dt <- as.data.table(detections)
    positions_dt$pos_ts <- positions_dt$up_timestamp_utc
    setkey(positions_dt, up_timestamp_utc)
    setkey(detections_dt, det_ts)
    detpos <- as_tibble(positions_dt[detections_dt, , roll = "nearest"])
    detpos_clean <- detpos %>%
      rename(dets_ts = up_timestamp_utc) %>%
      filter(abs(as.numeric(difftime(dets_ts, pos_ts, units = "secs"))) < PAR_DET_POS_TIMEDIFF) %>%
      mutate(dets_ts_5min = round_date(x = dets_ts, unit = "5 mins"))
      
    #rolljoin detpos with thermocline
    detpos_clean_dt <- as.data.table(detpos_clean)
    thermocline_wide_dt <- as.data.table(thermocline_wide)
    thermocline_wide_dt$interval_tmp <- thermocline_wide_dt$interval
    
    setkey(detpos_clean_dt, dets_ts)
    setkey(thermocline_wide_dt, interval_tmp)
    detpos_therm <- as_tibble(thermocline_wide_dt[detpos_clean_dt, roll = "nearest"])
    
    # clean and interpolate values
    #TODO: there is therm_part, account for that in interpolation, some (most) variables should use just center
    detpos_therm_interpolated <- detpos_therm %>% 
      rename(dets_ts = interval_tmp) %>% 
      filter(abs(as.numeric(difftime(dets_ts, interval, units = "secs"))) < 60*30) %>% #remove detections for which the thermocline log is further than 30 min
      rowwise() %>%
      mutate(det_therm_thickness = interpolate_thermocline_value_linear(logger_distances = c(distance_east, distance_west),
                                                                        logger_values = c(thickness_East, thickness_West),
                                                                        detection_distances = dist_from_zero),
             det_therm_temperature = interpolate_thermocline_value_linear(logger_distances = c(distance_east, distance_west),
                                                                          logger_values = c(temperature_East, temperature_West),
                                                                          detection_distances = dist_from_zero),
             det_therm_depth_center = interpolate_thermocline_value_linear(logger_distances = c(distance_east, distance_west),
                                                                           logger_values = c(depth_East, depth_West),
                                                                           detection_distances = dist_from_zero),
             det_therm_deviation = interpolate_thermocline_value_linear(logger_distances = c(distance_east, distance_west),
                                                                        logger_values = c(deviation_East, deviation_West),
                                                                        detection_distances = dist_from_zero),
             # TODO: det_therm_strengh
             )
      
  }
}