# Functions

#' Compute thermocline given depth and temperature vector
#' 
#' @param depth numeric vector of depths
#' @param temperature numeric vector of temperatures
#' @param diff_threshold threshold on slope (temperature/meter)
#' @param depth_res 
#' @return TODO:
#' @details depths are smoothed by monotonic function using `stats::splinefun()`. Thermocline is then detected on new curve.

compute_thermocline <- function(depth, temperature, diff_threshold = 2, depth_res = 0.1){
  #order the vectors in case it is not ordered
  if(length(depth) < 2) return(NULL)# do not return null - exception better
  if(!all(depth == cummax(depth))){
    depth <- depth[order(depth)]
    temperature <- temperature[order(depth)]
  }
  #create sequence of new depths
  depths_new <- seq(floor(min(depth)), ceiling(max(depth)), depth_res)
  #smooth the temperature profile
  frame_spline <- splinefun(x = depth, y = temperature, method = "monoH.FC", ties = mean)
  frame_new <- data.table(depth =  depths_new, temperature = frame_spline(depths_new))
  #compute derivations in each step (decrease of temperature per meter)
  frame_new[, temp_diff := c(diff(temperature)/depth_res, 0)]
  frame_new[, therm_diff_bool := -diff_threshold > temp_diff]
  #add also ending point from which the thermocline was not with diff.threshodl slope
  frame_new[shift(therm_diff_bool, n = 1, fill = F, type = "lag") == T & therm_diff_bool == F , therm_diff_bool := T]
  #add group when the validity of condition changed
  frame_new[, therm_split := rleid(therm_diff_bool)]
  #take only valid steps
  th_steps <- frame_new[therm_diff_bool == T]
  #assign new sequence of steps 1, 2, 3...
  th_steps[, step_order := rleid(therm_split)]
  th_steps[, step_order := abs(step_order-max(step_order))+1]
  #aggregate
  th_steps.agg <- th_steps[, .(depth_start = min(depth),
                               depth_end = max(depth),
                               temperature_start = temperature[1],
                               temperature_end = tail(temperature, 1),
                               depth_crit = depth[which(temp_diff == min(temp_diff))[1]],
                               temperature_crit = temperature[which(temp_diff == min(temp_diff))[1]]),
                           by = step_order]
  return(th_steps.agg)
}


#' Apply rolling function by time span
#' 
#' @param x vector of values
#' @param times sorted time vector - same length as x
#' @param span width of window in seconds 
roll_time_window <- function(x, times, span, FUN){
x_out <- vector(mode = class(x), length = length(x))
x_out[1:length(x_out)] <- NA
if(is.unsorted(times)){
  x <- x[order(times)]
  times <- times[order(times)]
}
FUN <- match.fun(FUN)
times_minus <- times - span
times_plus <- times + span
for(i in 1:length(x_out)){
  x_out[i] <- FUN(x[which(times > times_minus[i] & times < times_plus[i])])
}
return(x_out)
}
