# calculation of mean thermocline temperature for every single day (mean from both lines, calculated as diffrence between up_depth and down_depth)
thermocline <- here("data/products/thermocline_data.csv") %>%
  read_csv() %>%
  filter(step_order == 1 & slope == PAR_THERMOCLINE_SLOPE & therm_part == "center")

wind_data <- read_csv(here("data", "raw", "db", "wind_data.csv"))

wind_run <- wind_data %>% filter(parameter == "wind_run")

wind_direction <- wind_data %>% filter(parameter == "winddirection")

wind_speed <- wind_data %>% filter(parameter == "windspeed")


ggplot(wind_data, aes(x = md_timestamp_utc, y = val)) + 
  geom_point() +
  facet_wrap(~parameter, ncol = 1, scales = "free_y")

ggplot(wind_data, aes(x = md_timestamp_utc, y = val)) +
  geom_point() +
  facet_wrap(~parameter, ncol = 1, scales = "free_y")

ggplot(wind_direction, aes(x = val)) +
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

wind <- merge(wind_speed[, .(md_timestamp_utc, speed  = val)], wind_direction[, .(md_timestamp_utc, direction  = val)])
#azimuth of line is
#1.72522528309502 rad
#you can convert degrees into radian by
#(x/360)*2*pi 
wind[, direction.rad :=  (direction/360)*2*pi]

wind[, windproj := projectWindVector(angleOfLine = 1.72522528309502, distance = speed, direction = direction.rad)]
wind[, windproj.abs := abs(windproj)]
###SHINY PLOTTING ####

  time_range <- range(thermocline$interval)
  gp1 <- thermocline %>%
    filter(therm_part == "center" & step_order == 1) %>%
    ggplot(mapping = aes(x = interval , y = deviation, col = location))+
    geom_line() +
    guides(col = F) +
    xlab("Date") + 
    ylab("Deviation (m)")+  theme_minimal() + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  

  gp2 <- wind_speed %>%
    filter(md_timestamp_utc %between% time_range) %>%
    ggplot(aes(x = md_timestamp_utc, y = val)) + 
    geom_line() +
    xlab("Date") + ylab("Wind speed (m/s)")+  theme_minimal() + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
  gp3 <- wind_direction %>%
    filter(md_timestamp_utc %between% time_range) %>%
    ggplot(aes(x = md_timestamp_utc, y = val)) + 
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
