# Functions
# Functions -----------------------------------------------



#' Compute thermocline for profile using rLakeAnalyzer::wtr.layer function
#' @param x data.frame with two columns depth and temperature representing one time frame

compute_thermocline_wtr_layer <- function(x){
  wl <-  wtr.layer(depth = x$depth, measure = x$temperature, zmax = 25, thres = 0.1, z0 = 4)
  temperature_center <- approx(x = x$depth, y = x$temperature, xout = wl$cline)$y
  segments <- wl$segments[[1]]
  thermocline_segment_idx <- findInterval(vec = segments$segment_depth, x = wl$cline)
  
  tibble(
    depth_start = wl$mld,
    depth_crit = wl$cline,
    depth_end = segments[thermocline_segment_idx + 1,]$segment_depth,
    temperature_start = segments[thermocline_segment_idx,]$segment_measure,
    temperature_crit = temperature_center,
    temperature_end = segments[thermocline_segment_idx + 1,]$segment_measure
  )  %>%
    mutate(thickness = depth_end - depth_start, 
           strength = temperature_start - temperature_end) %>%
    mutate(mean_gradient = strength / thickness)
}



#' Compute thermocline for profile using rLakeAnalyzer::thermo.depth and rLakeAnalyzer::meta.depths functions
#' @param x data.frame with two columns depth and temperature representing one time frame

compute_thermocline_thermo_depth <- function(x){
  slope <- 0.2
  seasonal <- T
  
  depth_crit <- thermo.depth(wtr = x$temperature, x$depth, seasonal = seasonal)
  temperature_crit <- approx(x = x$depth,
                               y = x$temperature,
                               xout = depth_crit)$y
  
  meta_depths <- meta.depths(wtr = x$temperature,
                             x$depth, slope = slope,
                             seasonal = seasonal)
  meta_temperatures <- approx(x = x$depth,
                              y = x$temperature,
                              xout = meta_depths)$y
  
  tibble(
    depth_start = meta_depths[1],
    depth_crit = depth_crit,
    depth_end = meta_depths[2], 
    temperature_start = meta_temperatures[1],
    temperature_crit = temperature_crit,
    temperature_end = meta_temperatures[2]) %>%
    mutate(thickness = depth_end - depth_start, 
           strength = temperature_start - temperature_end) %>%
    mutate(mean_gradient = strength / thickness)
}

#' Compute thermocline given depth and temperature vector
#'
#' @param depth numeric vector of depths
#' @param temperature numeric vector of temperatures
#' @param diff_threshold threshold on slope (temperature/meter)
#' @param depth_res
#' @return TODO:
#' @details depths are smoothed by monotonic function using `stats::splinefun()`. Thermocline is then detected on new curve.

compute_thermocline <- function(depth, temperature, diff_threshold = 2, depth_res = 0.1) {
  # order the vectors in case it is not ordered
  if (length(depth) < 3) stop("Cannot compute thermocline with less than 3 points")
  if (!all(depth == cummax(depth))) {
    depth <- depth[order(depth)]
    temperature <- temperature[order(depth)]
  }
  # create sequence of new depths
  frame_new <- as.data.table(interpolate_temperature_profile(depth, temperature))
  # compute derivations in each step (decrease of temperature per meter)
  frame_new[, therm_diff_bool := -diff_threshold > slope]
  # add also ending point from which the thermocline was not with diff.threshodl slope
  frame_new[data.table::shift(therm_diff_bool, n = 1, fill = F, type = "lag") == T & therm_diff_bool == F, therm_diff_bool := T]
  # add group when the validity of condition changed
  frame_new[, therm_split := rleid(therm_diff_bool)]
  # take only valid steps
  th_steps <- frame_new[therm_diff_bool == T]
  if (nrow(th_steps) == 0) {
    th_steps.agg <- th_steps[, .(
      step_order = 1,
      depth_start = as.numeric(NA),
      depth_end = as.numeric(NA),
      temperature_start = as.numeric(NA),
      temperature_end = as.numeric(NA),
      depth_crit = as.numeric(NA),
      temperature_crit = as.numeric(NA),
      thickness = NA,
      strength = NA,
      mean_gradient = NA
    )]
  } else {
    # assign new sequence of steps 1, 2, 3...
    th_steps[, step_order := rleid(therm_split)]
    th_steps[, step_order := abs(step_order - max(step_order)) + 1]
    # aggregate
    th_steps.agg <- th_steps[, .(
      depth_start = min(depth),
      depth_end = max(depth),
      temperature_start = temperature[1],
      temperature_end = tail(temperature, 1),
      depth_crit = depth[which(slope == min(slope))[1]],
      temperature_crit = temperature[which(slope == min(slope))[1]],
      thickness = max(depth) - min(depth),
      strength = temperature[1] - tail(temperature, 1),
      mean_gradient =  (temperature[1] - tail(temperature, 1)) / (max(depth) - min(depth))
    ),
    by = step_order
    ]
  }
  return(th_steps.agg)
}

compute_thermocline_dplyr <- function(depth, temperature, diff_threshold = 2, depth_res = 0.1) {
  # order the vectors in case it is not ordered
  if (length(depth) < 3) stop("Cannot compute thermocline with less than 3 points")
  if (!all(depth == cummax(depth))) {
    depth <- depth[order(depth)]
    temperature <- temperature[order(depth)]
  }
  # create sequence of new depths
  depths_new <- seq(floor(min(depth)), ceiling(max(depth)), depth_res)
  # smooth the temperature profile
  temperature_model <- splinefun(x = depth, y = temperature, method = "monoH.FC", ties = mean)
  temperatures_new <- temperature_model(depths_new)
  # compute derivations in each step (decrease of temperature per meter)
  temperature_diff <- c(diff(temperatures_new) / depth_res, 0)
  is_thermocline <- -diff_threshold > temperature_diff
  # add also ending point from which the thermocline was not with diff.threshodl slope
  is_thermocline[data.table::shift(is_thermocline, n = 1, fill = F, type = "lag") == T] <- T

  if (!any(is_thermocline)) {
    result <- tibble(
      step_order = 1,
      depth_start = as.numeric(NA),
      depth_end = as.numeric(NA),
      temperature_start = as.numeric(NA),
      temperature_end = as.numeric(NA),
      depth_crit = as.numeric(NA),
      temperature_crit = as.numeric(NA)
    )
  } else {
    # add group when the validity of condition changed
    thermocline_steps <- rleid(is_thermocline)
    # Start with N. 1 at from bottom up

    thermoclines_tb <- dplyr::tibble(
      step_order = rleid(thermocline_steps[is_thermocline]),
      temperature = temperatures_new[is_thermocline],
      depth = depths_new[is_thermocline],
      temperature_diff = temperature_diff[is_thermocline]
    ) %>% mutate(step_order = (max(step_order) - step_order) + 1)

    # aggregate
    result <- thermoclines_tb %>%
      group_by(step_order) %>%
      summarize(
        depth_start = min(depth),
        depth_end = max(depth),
        temperature_start = temperature[depth == min(depth)][1],
        temperature_end = temperature[depth == max(depth)][1],
        depth_crit = depth[temperature_diff == min(temperature_diff)][1],
        temperature_crit = temperature[temperature_diff == min(temperature_diff)][1],
        .groups = "keep"
      )
  }
  return(as.data.table(result))
}

#' Apply rolling function by time span
#'
#' @param x vector of values
#' @param times sorted time vector - same length as x
#' @param span width of window in seconds
roll_time_window <- function(x, times, span, FUN) {
  x_out <- vector(mode = class(x), length = length(x))
  x_out[1:length(x_out)] <- NA
  if (is.unsorted(times)) {
    x <- x[order(times)]
    times <- times[order(times)]
  }
  FUN <- match.fun(FUN)
  times_minus <- times - span
  times_plus <- times + span
  for (i in 1:length(x_out)) {
    x_out[i] <- FUN(x[which(times > times_minus[i] & times < times_plus[i])])
  }
  return(x_out)
}

#' Functions loading data
load_hobo_data <- function() {
  read_csv(here("data", "raw", "db", "hobo_data.csv"), col_types = c("dTdddccdTTTd"))
}

load_temperature_data <- function() {
  read_csv(here("data", "products", "temperature_data.csv"), col_types = c("cTddd"))
}

load_logger_positions <- function() {
  here("data/raw/db/logger_positions.shp") %>% 
    st_read(quiet = T) %>%
    rename(location = locatin,
           dist_from_zero = dst_fr_,
           location_order = lctn_rd)
}

load_wind_data <- function() {
  read_csv(here("data", "raw", "db", "wind_data.csv"), col_types = c("cccccddTd"))
}

#' Smooth temperature profile
#'
#' @inheritParams compute_thermocline
#' @details
interpolate_temperature_profile <- function(depth, temperature, depth_res = 0.1) {
  if (!all(depth == cummax(depth))) {
    stop("Depth must be in decreasing order")
  }
  depths_new <- seq(floor(min(depth)), ceiling(max(depth)), depth_res)
  temperature_monotonic <- cummin(temperature)
  temperature_profile_fun <- splinefun(x = depth, y = temperature_monotonic, method = "hyman", ties = mean)
  temperature_new <- temperature_profile_fun(depths_new)
  return(list(depth = depths_new, temperature = temperature_new, slope = c(diff(temperature_new) / depth_res, 0)))
}



#' Get day/night for given time
#'
#' @param x vector of times
#' @param lat latitude
#' @param lon longitude
#' @param label_day character vector of length 1 - label for day
#' @param label_night character vector of length 1 - label for day
#'
#' @return character vector giving night or day
#' @export
#'
#' @examples
#' get_diel_period(x = seq.POSIXt(from = Sys.time(), to = Sys.time() + 84600, length.out = 20))
get_diel_period <- function(x, lat = 49.5765639, lon = 14.6637706, label_day = "day", label_night = "night") {
  if (length(x) == 0) stop("Length of the input is 0")
  if (all(is.na(x))) {
    return(as.character(x))
  }
  sunset_sunrise <- get_sunsets_sunrises(x, lat = lat, lon = lon)
  day_nigth <- as.character(ifelse(x > sunset_sunrise$sunrise & x < sunset_sunrise$sunset, label_day, label_night))
  return(day_nigth)
}


#' Get sunset and sunrise times for each given time
#'
#' For each given time, function returns sunset and sunrise time on that given day date
#'
#' @inheritParams get_diel_period
get_sunsets_sunrises <- function(x, lat = 49.5765639, lon = 14.6637706) {
  sunrise_sunset <- get_sunset_sunrise(x, lat, lon)
  sunrise_sunset$Date <- as.Date(sunrise_sunset$sunrise)
  res <- merge(data.frame(x, xorder = 1:length(x), Date = as.Date(x)), sunrise_sunset, by = "Date", all.x = T)
  res[order(res$xorder), c("sunset", "sunrise")]
}

#' Get sunset and sunrise times for time span given by time vector
#'
#' for time span given by time vector, function returns sunset and sunrise times
#'
#' @inheritParams get_diel_period
#' @importFrom StreamMetabolism sunrise.set
get_sunset_sunrise <- function(x, lat = 49.5765639, lon = 14.6637706) {
  from <- min(x, na.rm = T)
  to <- max(x, na.rm = T)
  StreamMetabolism::sunrise.set(lat, lon, from, num.days = difftime(time1 = to, time2 = from, units = "days") + 2)
}

#' Get night time polygons
#'
#' @inheritParams get_diel_period
#'
#' @return
#' @export
#'
#' @examples
get_nighttime_polygons <- function(x, lat = 49.5765639, lon = 14.6637706) {
  sunsetsunrise <- get_sunset_sunrise(x, lat = lat, lon = lon)
  nightime.vector <- sort(c(sunsetsunrise[-1, 1], sunsetsunrise[-nrow(sunsetsunrise), 2]))
  x <- rep(nightime.vector, each = 2)
  y <- rep(c(-Inf, Inf), times = length(nightime.vector))
  xord <- rep(c(1, 2, 4, 3), times = length(nightime.vector) / 2)
  night_polygons <- data.frame(x,
    y,
    id = rep(1:(round(length(nightime.vector) / 2)), each = 4),
    xord = xord
  )
  night_polygons <- night_polygons[order(xord), ]
  return(night_polygons)
}


#' For given positions, get distance of their projections from point zero (one end of lake axis/shore)
#'
#' @param positions spatial points of positions
#' @param lake_axis spatial line to project to
#' @param point_zero spatial point on one side of the lake to measure the distance from
#' @return numeric vector of distances for each given position's projection on lake_axis
compute_distance_from_point_zero <- function(positions, lake_axis, point_zero) {
  nearest_points <- st_nearest_points(positions, lake_axis) %>% st_cast("POINT")
  # previous step returned both starts and ends. Pick only points on axis
  projected_positions <- nearest_points[seq(2, length(nearest_points), 2)]
  as.numeric(st_distance(projected_positions, point_zero))
}



#' Intepolation of thermocline value in an arbitrary distance.
#' @param logger_distances vector of logger distances form point zero
#' @param logger_values vector of values measured on loggers
#' @param detection_distances vector of detection distances from point zero
#' @details If the length of logger_values and detection_distances is equal, the values are assumed to be matching.
#' @return vector of values at place of detections
interpolate_thermocline_value_linear <- function(logger_distances, logger_values, detection_distances) {
  approx(x = logger_distances, y = logger_values, xout = detection_distances, rule = 2)$y
}


#' Function joining positions to each detection
#' @param detections df of detections
#' @param positions df of positions
#' @details
join_detections_positions <- function(detections, poisitions, max_timediff) {
  # rolljoin positions and detections
  positions <- data.table::as.data.table(positions)
  detections <- data.table::as.data.table(detections)
  positions$pos_ts <- positions$up_timestamp_utc
  data.table::setkey(positions, up_timestamp_utc)
  data.table::setkey(detections, det_ts)
  detpos <- tidyr::as_tibble(positions[detections, , roll = "nearest"])
  detpos_clean <- detpos %>%
    rename(dets_ts = up_timestamp_utc) %>%
    filter(abs(as.numeric(difftime(dets_ts, pos_ts, units = "secs"))) < max_timediff) %>%
    mutate(dets_ts_5min = round_date(x = dets_ts, unit = "5 mins"))
}

#' Function joining detection records to thermocline
#' @param detections df of detections (dets_ts, ...)
#' @param thermocline df of thermocline (therm_part, thermocline_ts, ...)
#' @details for each thermocline part, the detections are joined separately -> result can contain more records of one detections
#' @return data.frame of joined detections to thermocline
join_detections_thermocline <- function(detections, thermocline, max_timediff) {
  x <- merge(as.data.frame(detections), data.frame(therm_part = unique(thermocline$therm_part)), all = T)
  x_dt <- as.data.table(x)
  # x_dt$thermocline_ts <- x_dt$dets_ts
  y_dt <- as.data.table(thermocline)
  y_dt$dets_ts_5min <- y_dt$thermocline_ts
  # y_dt$interval_tmp <- y_dt$thermocline_ts
  setkey(x_dt, therm_part, dets_ts_5min)
  setkey(y_dt, therm_part, dets_ts_5min)
  detpos_therm <- as_tibble(y_dt[x_dt, roll = "nearest"]) %>%
    filter(abs(as.numeric(difftime(dets_ts_5min, thermocline_ts, units = "secs"))) < max_timediff) # remove detections for which the thermocline log is further than 30 min
}
