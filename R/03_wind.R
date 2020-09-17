# calculation of mean thermocline temperature for every single day (mean from both lines, calculated as diffrence between up_depth and down_depth)
thermocline <- thermocline_location


ggplot(wind.long, aes(x = md_timestamp_utc, y = val)) + 
  geom_point() +
  facet_wrap(~parameter, ncol = 1, scales = "free_y")

ggplot(wind.long, aes(x = md_timestamp_utc, y = val)) +
  geom_point() +
  facet_wrap(~parameter, ncol = 1, scales = "free_y")

ggplot(winddirection, aes(x = val)) +
  geom_histogram(binwidth = 5) +
  coord_polar() 



#function projecting one vector to a line (vector version)
#all angles has to be positive in radians!
#angleOfLine - azimuth of logger line
#distance - measure of distance (it can be windspeed, windrun or whatever)
projectWindVector <- function(angleOfLine, distance, direction){
  wind.vectors <- cbind(distance * cos(distance), distance * sin(distance))
  wind.projected <- wind.vectors * cos(direction - angleOfLine)
  res <- sign(cos((direction-angleOfLine))) * sqrt(wind.projected[,1]^2 + wind.projected[,2]^2)
  return(res)
}

wind <- merge(windspeed[, .(md_timestamp_utc, speed  = val)], winddirection[, .(md_timestamp_utc, direction  = val)])
#azimuth of line is
#1.72522528309502 rad
#you can convert degrees into radian by
#(x/360)*2*pi 
wind[, direction.rad :=  (direction/360)*2*pi]

wind[, windproj := projectWindVector(angleOfLine = 1.72522528309502, distance = speed, direction = direction.rad)]
wind[, windproj.abs := abs(windproj)]
###SHINY PLOTTING ####

  time_range <- range(thermocline$interval)
  gp1 <- ggplot(data = th_m[therm_part == "center" & step_order == 1 ], 
                mapping = aes(x = interval , y = deviation , col = location))+
    geom_line() +
    guides(col = F) +
    xlab("Date") + 
    ylab("Deviation (m)")+  theme_minimal() + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  

  gp2 <- ggplot(data = wind[md_timestamp_utc %between% time_range], aes(x = md_timestamp_utc, y = speed)) + 
    geom_line() +
    xlab("Date") + ylab("Wind speed (m/s)")+  theme_minimal() + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
  gp3 <- ggplot(data = wind[md_timestamp_utc %between% time_range], aes(x = md_timestamp_utc, y = direction)) + 
    geom_line() +
    scale_y_continuous( expand = c(0, 0), breaks = c(0, 90, 180, 270, 360), labels = c("N", "E", "S", "W", "N"), limits = c(0,360)) +
    xlab("Date") + ylab("Wind direction")  +  theme_minimal() 
  g1 <- ggplotGrob(gp1)
  g2 <- ggplotGrob(gp2)
  g3 <- ggplotGrob(gp3)
  grid.newpage()
  grid.draw(rbind(g1, g2, g3, size = "last"))
  
# Data export --------------------------------------
# TODO:  wind strength (mean, median, range - split into day and night)
