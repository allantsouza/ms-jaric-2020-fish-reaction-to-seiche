############################################### ###
# Setting of dataset and libraries ####
#setwd("~/Macfish/Seche")
#load packages
library(RODBC)
library(lubridate)
library(StreamMetabolism) # for calculation of sunset and sunrise 
library(rgdal)
library(SDMTools)
library(rgeos)
library(maptools)
library(lubridate)
library(sp)    # spatial data 
library(doMC)  # this package and next one are needed for processing data on more cores
library(parallel)
library(geosphere)
library(ggplot2)
library(adehabitatHR)
library (ggthemes) # extension of ggplots, used for home range visualsation
library(leaflet) # visualization of home range on the google maps, likely not needed
library (lunar) # needed for lunar phase calculation
library(data.table)
library(postGIStools)
library(RPostgreSQL)
library(dygraphs)


thermpart_par <- "center"
step_order_par <- 1
degree_per_meter_par <- 1

DateRange <- range("2015-06-16","2015-10-15")


######################################
con <-  dbConnect(drv = PostgreSQL(), dbname ="macfishdb", host="localhost", user= "macfishuser", password = "macf1shus3r!")


#load thermocline dataset
select.templine <-paste("
                        SELECT a.th_timestamp_utc as timestamp_utc, th_step_order step_order, th_degree_per_meter degree_per_meter, th_thermpart thermpart, 
                        th_depth, th_temperature temperature, a.lp_pos_name pos_name, b.lt_lake lake, b.lp_pos_order pos_order FROM at_macfish.thermocline a INNER JOIN at_macfish.loggerpos b ON a.lp_pos_name = b.lp_pos_name 
                        AND lt_lake = 'Chabarovice';", sep = "") 

# calculation of mean thermocline temperature for every single day (mean from both lines, calculated as diffrence between up_depth and down_depth)
thermocline <- data.table(dbGetQuery(con, select.templine, stringsAsFactors = F))
thermocline[, Date := as.Date(timestamp_utc)]

#DateRange <- range(thermocline$date)
thermocline.sub <- thermocline[Date <= DateRange[2] & Date >= DateRange[1] & thermpart == thermpart_par & degree_per_meter == degree_per_meter_par & step_order == step_order_par ]
ggplot(thermocline.sub , aes(x = timestamp_utc, y = th_depth, col = pos_name))+geom_line()+geom_smooth()

#Tilt of thermocline calculation (creation of "shift.dataset" datatable)
setkey(thermocline.sub, pos_order, timestamp_utc)
therm.mean <- thermocline.sub[,.(mean.depth = mean(th_depth), d.count = .N),by = .(lake, timestamp_utc,step_order, degree_per_meter,thermpart)]

#therm.mean[, mean.depth.smooth := rollapply(width = 12*24, FUN = mean, data = mean.depth, fill = NA), by = .( lake, step_order, degree_per_meter, thermpart)]
setkey(therm.mean,  timestamp_utc)
therm.mean[, mean.depth.smooth := rollTimeWindow(win.width = 60*60*24*4, FUN = mean, x = mean.depth, time.vec = timestamp_utc), by = .( lake, step_order, degree_per_meter, thermpart)]

#therm mean now contains depth of thermocline in balanced state without wind
therm.mean

#load fish detections
select.fish.qu <- 
  "SELECT tf_tag_sn
FROM at_macfish.taggedfish
WHERE
tf_species = 'tench' AND tf_lake = 'Chabarovice'
;
"
# picking appropriate tag ids, needed for all types of query above
tag_ids <- dbGetQuery(con, select.fish.qu)
tag_ids <- tag_ids[,1]
#get tagids 
select.detections <-paste("
                          WITH posth AS (SELECT th.*, lake, up_tag_sn, up_disttoshore, pos.up_timestamp_utc, pos.up_depth FROM 
                          (SELECT up_tag_sn, up_id, up_disttoshore, up_timestamp_utc, up_depth, up_lake lake FROM at_macfish.umap_pos_powfilter WHERE up_validpos AND fishvalid AND up_tag_sn IN ('", paste(tag_ids, collapse = "','", sep = ""),"')) pos INNER JOIN at_macfish.umap_pos_thermocline th ON pos.up_id = th.up_id AND pth_thermpart = '", thermpart_par,"')
                          
                          SELECT pos.*, dd_depth, dd_timestamp_utc FROM at_macfish.detsdepth dd 
                          INNER JOIN (SELECT posth.*, b.dd_id FROM posth INNER JOIN at_macfish.detsdepth_to_umap_pos b ON b.up_id = posth.up_id) pos ON
                          pos.dd_id = dd.dd_id;", sep = "")

detections <- data.table(dbGetQuery(con, select.detections, stringsAsFactors = F))

detections[, det_pos_timediff := abs(as.numeric(difftime(dd_timestamp_utc, up_timestamp_utc, units = "secs")))]
#remove detections which are further than 30 min from closest position
detections <- detections[det_pos_timediff < 60*15]


#load pos_thermocline and add site-specific thermocline depth to each detection - roll join
detections[, timestamp_utc := dd_timestamp_utc]
detections[, thermpart := pth_thermpart]
detections[, step_order := pth_step_order]
detections[, degree_per_meter := pth_degree_per_meter]
therm.mean[, th_timestamp_utc := timestamp_utc]
#compute day and night
detections[, day_night :=  getNightDay(timestamp_utc)]
detections[, timestamp5m := round.POSIXct(dd_timestamp_utc, 5*60)]

#rolljoin
setkey(detections, lake, step_order, degree_per_meter, thermpart, timestamp_utc)
setkey(therm.mean, lake, step_order, degree_per_meter, thermpart, timestamp_utc)
detections.m <- therm.mean[detections, roll = "nearest"]
#remove detections for which the thermocline log is further than 30 min
detections.m[!is.na(th_timestamp_utc), therm_det_timediff := abs(as.numeric(difftime(dd_timestamp_utc, th_timestamp_utc, units = "secs")))]
detections.m.sub <- detections.m[!is.na(th_timestamp_utc) & therm_det_timediff < 60*15]
#make difference between fish depth and site-specific thermocline depth
detections.m.sub[, fish_minus_therm_depth := dd_depth - pth_thermocline_depth]
#compute how much is thermocline tilted - out of balance
detections.m.sub[, therm_minus_balance :=  pth_thermocline_depth - mean.depth.smooth]
detections.m.sub[, fish_minus_balance :=  dd_depth - mean.depth.smooth]


#compare situations when the thermocline is up and when down
therm_minus_balance_lower <- -1
therm_minus_balance_upper <- 1
detections.m.sub[, therm_minus_balance_cat :=  factor(findInterval(x = therm_minus_balance, vec = c(therm_minus_balance_lower, therm_minus_balance_upper)),
                                                      levels = 0:2, labels = c(paste("<", therm_minus_balance_lower), paste("no-disbalance"), paste(">", therm_minus_balance_upper)))]

# ggplot(detections.m.sub[], aes(x = therm_minus_balance_cat, y = fish_minus_balance, group = interaction(therm_minus_balance_cat, day_night), col = day_night))+geom_boxplot()+facet_wrap(~up_tag_sn, scales = "free")

#DATA TABLE NAME = detections.m.sub
#NAMES OF COLUMNS FOR KEY PARAMETERS:
# fish depth = dd_depth
# disbalanced = therm_minus_balance
# FishID = up_tag_sn
# Diel = day_night

#LATER ON, IN CASE WE NEED MEAN DEPTH OF A THERMOCLINE (WITHOUT SEICHE) AS AN INDICATION OF SEASONS:
# thermocline depth = mean_depth_smooth



