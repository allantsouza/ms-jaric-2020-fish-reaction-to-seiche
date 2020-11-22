# Script for processing raw detections and positions and merging with thermocline
# Assigning rough lateral position (distance from Western shore) to each detection and merging with thermocline. 
# OUTPUT: detections with distance from shore and interpolated values of thermocline in their specific location


# Loading data ------------------------------------------------------------


fish_raw <- read_csv(file = "data/raw/fishIDs.csv", col_types = "ccdc")
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
  read_csv(col_types = c("thermocline_ts" = "T"))

thermocline_wide <- thermocline %>% 
  pivot_wider(id_cols = c("thermocline_ts", 
                          "therm_part",
                          "lake_therm_depth_smoothed",
                          "lake_therm_temperature_smoothed"), names_from = "location",
              values_from = c("deviation", "thickness", "temperature", "depth", "strength")) %>%
  group_by(thermocline_ts) %>%
  filter(!is.na(depth_East) & !is.na(depth_West)) %>%
  ungroup()
  




# For each fish -----------------------------------------------------------


for(i in 1:length(tag_sns)){
  # load positions and detections - skip if the file does not exist
  position_i_fp <- here(position_dir, paste0(tag_sns[i], ".csv"))
  if(!file.exists(position_i_fp)){
    next
  }
  positions <- read_csv(file = position_i_fp, col_types = c(up_depth = "d")) %>%
    st_as_sf(coords = c("x", "y"), crs = 32633)
  
  
  detections_i_fp <- here(detection_dir, paste0(tag_sns[i], ".csv"))
  if(!file.exists(position_i_fp)){
    next
  }
  detections <- read_csv(file = detections_i_fp, col_types = c("dcdT"))
  
  # get distance from shore
  positions$dist_from_zero <- compute_distance_from_point_zero(positions = positions,
                                                               lake_axis = lake_axis,
                                                               point_zero = point_zero)
    
  # rolljoin detections to positions
  detpos <- join_detections_positions(detections, positions, max_timediff = PAR_DET_POS_TIMEDIFF)

  # rolljoin detpos with thermocline (do not join more than 30 mins appart)
  detpos_therm <- join_detections_thermocline(detpos, thermocline_wide, max_timediff = 60*30)
  if(nrow(detpos_therm) == 0){
    next
  }
  
  # clean and interpolate values into place of detection
  detpos_therm_interpolated <- detpos_therm %>% 
    rowwise() %>% # interpo
    mutate(det_therm_thickness = interpolate_thermocline_value_linear(logger_distances = c(distance_east, distance_west),
                                                                      logger_values = c(thickness_East, thickness_West),
                                                                      detection_distances = dist_from_zero),
           det_therm_temperature = interpolate_thermocline_value_linear(logger_distances = c(distance_east, distance_west),
                                                                        logger_values = c(temperature_East, temperature_West),
                                                                        detection_distances = dist_from_zero),
           det_therm_depth = interpolate_thermocline_value_linear(logger_distances = c(distance_east, distance_west),
                                                                         logger_values = c(depth_East, depth_West),
                                                                         detection_distances = dist_from_zero),
           det_therm_gradient = interpolate_thermocline_value_linear(logger_distances = c(distance_east, distance_west),
                                                                  logger_values = c(depth_East, depth_West),
                                                                  detection_distances = dist_from_zero),
           ) %>%
    mutate(det_therm_deviation = lake_therm_depth_smoothed - det_therm_depth)
  
  # Widen so all therm_parts in separate columns
  detpos_therm_interpolated_wide <- detpos_therm_interpolated %>%
    pivot_wider(id_cols = c("tag_sn", "dets_ts", "det_depth", "thermocline_ts", "dist_from_zero"), 
              names_from = "therm_part", 
              values_from = c("lake_therm_depth_smoothed", "det_therm_temperature", "det_therm_depth", "det_therm_deviation", "depth_East", "depth_West")) %>%
    mutate(det_therm_strength = det_therm_temperature_start - det_therm_temperature_end) %>%
    mutate(lake_disbalance = depth_East_center - depth_West_center) %>%
    mutate(is_valid_seiche = abs(lake_disbalance) > 0.69999) %>%
    mutate(lake_therm_thickness_smoothed = lake_therm_depth_smoothed_end - lake_therm_depth_smoothed_start) %>%
    mutate(diel_period = get_diel_period(dets_ts))

  # Export
  write_csv(x = detpos_therm_interpolated_wide, path = here("data/products/fish/", paste0(tag_sns[i], ".csv")))
}
