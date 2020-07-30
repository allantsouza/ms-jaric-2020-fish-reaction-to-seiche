

#Data loading script
con <-  dbConnect(drv = PostgreSQL(), dbname ="macfishdb", host="172.21.3.20", user= "macfishuser", password = "macf1shus3r!")

#load thermocline dataset
select.templine <-paste("
                        SELECT a.th_timestamp_utc as timestamp_utc, th_step_order step_order, th_degree_per_meter degree_per_meter, th_thermpart thermpart, 
                        th_depth, th_temperature temperature, a.lp_pos_name pos_name, b.lt_lake lake, b.lp_pos_order pos_order FROM at_macfish.thermocline a INNER JOIN at_macfish.loggerpos b ON a.lp_pos_name = b.lp_pos_name 
                        AND lt_lake = 'Chabarovice';", sep = "") 

# calculation of mean thermocline temperature for every single day (mean from both lines, calculated as diffrence between up_depth and down_depth)
thermocline <- data.table(dbGetQuery(con, select.templine, stringsAsFactors = F))
thermocline[, Date := as.Date(timestamp_utc)]
thermocline.sub <- thermocline[Date <= DateRange[2] & Date >= DateRange[1]  & degree_per_meter == degree_per_meter_par & step_order == step_order_par ]
# ggplot(thermocline.sub[], aes(x = timestamp_utc, y = th_depth , col = thermpart))+
#   geom_line()+facet_wrap(~pos_name, ncol = 1)
# 

#Tilt of thermocline calculation (creation of "shift.dataset" datatable)
  setkey(thermocline.sub, pos_order, timestamp_utc)
  #get mean depth of thremocline in each 5min interval
  therm.mean <- thermocline.sub[,.(mean.therm.depth = mean(th_depth)), by = .(lake, timestamp_utc,step_order, degree_per_meter, thermpart)]
  #smooth mean depth by 4 days moving window
  setkey(therm.mean,  timestamp_utc)
  therm.mean[, balanced_therm_depth := rollTimeWindow(win.width = therm_bal_ws, FUN = mean, x = mean.therm.depth, time.vec = timestamp_utc), by = .(lake, step_order, degree_per_meter, thermpart)]
  therm.mean[, "mean.therm.depth" := NULL]
#therm mean now contains depth of thermocline in balanced state without wind

  # ggplot(therm.mean[], aes(x = timestamp_utc, y = balanced_therm_depth, col = thermpart))+
  #   geom_line()
#load fish detections
select.fish.qu <- paste("SELECT tf_tag_sn FROM at_macfish.taggedfish 
WHERE tf_species = '", species,"' AND tf_lake = 'Chabarovice'", sep = "")
# picking appropriate tag ids, needed for all types of query above
tag_ids <- dbGetQuery(con, select.fish.qu)
tag_ids <- tag_ids[,1]
#get tagids 
select.detections <-paste("
                          WITH posth AS (SELECT th.*, lake, up_tag_sn, up_disttoshore, up_timestamp_utc FROM 
                          (SELECT up_tag_sn, up_id, up_disttoshore, up_timestamp_utc, up_depth, up_lake lake FROM at_macfish.umap_pos_powfilter WHERE
                              up_validpos AND fishvalid AND up_tag_sn IN ('", paste(tag_ids, collapse = "','", sep = ""),"')) pos 
                          INNER JOIN 
                            (SELECT * FROM at_macfish.umap_pos_thermocline WHERE pth_step_order = ", step_order_par, " AND pth_degree_per_meter = ", degree_per_meter_par, " 
                                                                                AND pth_thermpart IN ('", paste(c("center","start","end"), collapse = "','"),"')
                              ) as th
                          ON pos.up_id = th.up_id )
                          
                          SELECT pos.*, dd_depth, dd_timestamp_utc FROM at_macfish.detsdepth dd 
                          INNER JOIN (SELECT posth.*, b.dd_id FROM posth INNER JOIN at_macfish.detsdepth_to_umap_pos b ON b.up_id = posth.up_id) pos ON
                          pos.dd_id = dd.dd_id;", sep = "")
detections <- data.table(dbGetQuery(con, select.detections, stringsAsFactors = F))

detections[, det_pos_timediff := abs(as.numeric(difftime(dd_timestamp_utc, up_timestamp_utc, units = "secs")))]
#remove detections which are further than 30 min from closest position
detections <- detections[det_pos_timediff < excl_dets_further_than_pos]


#load pos_thermocline and add site-specific thermocline depth to each detection - roll join
detections[, timestamp_utc := dd_timestamp_utc]
detections[, thermpart := pth_thermpart]
detections[, step_order := pth_step_order]
detections[, degree_per_meter := pth_degree_per_meter]
detections[,Date := as.Date(dd_timestamp_utc)]
therm.mean[, th_timestamp_utc := timestamp_utc]

detections[, dd_timestamp_5_min := round_POSIXct(dd_timestamp_utc, 5*60)]

#selecting only fish present for more than 30 days
no.days.available <- detections[, .(no.days = length(unique(Date))), by = up_tag_sn]
detections <- detections[up_tag_sn %in% no.days.available[no.days>30]$up_tag_sn]

#rolljoin
setkey(detections, lake, step_order, degree_per_meter, thermpart, timestamp_utc)
setkey(therm.mean, lake, step_order, degree_per_meter, thermpart, timestamp_utc)
detections.m <- therm.mean[detections, roll = "nearest"]
#remove detections for which the thermocline log is further than 30 min
detections.m[!is.na(th_timestamp_utc), therm_det_timediff := abs(as.numeric(difftime(dd_timestamp_utc, th_timestamp_utc, units = "secs")))]
detections.m.sub <- detections.m[!is.na(th_timestamp_utc) & therm_det_timediff < excl_dets_further_than_therm]

#make difference between fish depth and site-specific thermocline depth
#detections.m.sub[, fish_minus_therm_depth := dd_depth - pth_thermocline_depth]
#compute how much is thermocline tilted - out of balance
#detections.m.sub[, therm_minus_balance :=  pth_thermocline_depth - balanced_therm_depth]
#detections.m.sub[, fish_minus_balance :=  dd_depth - balanced_therm_depth]
oldw <- getOption("warn")
options(warn = -1)
detections.m.sub[, c("therm_det_timediff", "th_timestamp_utc", "up_id","pth_distance_from_loggerpos1", "up_timestamp_utc",
                     "dd_timestamp_5_min", "det_pos_timediff", "pth_step_order", "pth_degree_per_meter","step_order", "degree_per_meter", "timestamp_utc", 'pth_thermpart') := NULL]
setnames(detections.m.sub, c("pth_thermocline_temperature",  "pth_thermocline_depth"), c("thermocline_temperature", "thermocline_depth"))
options(warn = oldw)

detections.wide <- dcast(detections.m.sub, formula = lake +up_tag_sn + dd_depth + dd_timestamp_utc + dd_id+up_disttoshore ~ thermpart, value.var = c("balanced_therm_depth", "thermocline_temperature", "thermocline_depth"))

#compute day and night
detections.wide[, day_night :=  getNightDay(dd_timestamp_utc)]

#compare situations when the thermocline is up and when down
# therm_minus_balance_lower <- -1
# therm_minus_balance_upper <- 1
# detections.m.sub[, therm_minus_balance_cat :=  factor(findInterval(x = therm_minus_balance, vec = c(therm_minus_balance_lower, therm_minus_balance_upper)),
#                                                       levels = 0:2, labels = c(paste("<", therm_minus_balance_lower), paste("balanced"), paste(">", therm_minus_balance_upper)))]
# 
# detections.m.sub[, seiche_state :=  factor(findInterval(x = therm_minus_balance, 
#                                                         vec = c(therm_minus_balance_lower, therm_minus_balance_upper)),
#                                                       levels = 0:2, labels = c('downward', 'balanced', 'upward'))]



#Changing the name of data.table to apply simplifications and standarization of column names
##Modifications after data went to Souza
detections_std <- copy(detections.wide) 
setnames(detections_std, c("dd_timestamp_utc", "up_tag_sn", "up_disttoshore", "dd_depth"), 
         c("timestamp_utc", "fishid", "disttoshore", "fishdepth"))


fishids <- unique(detections_std$fishid)
fishid_shiny_def <- fishids[1]
time_range_shiny <- range(detections_std$timestamp_utc)
from.ss <- time_range_shiny[1]
to.ss <- time_range_shiny[2]
sunrise.sunset.t <- sunrise.set(lat = 50.653388, long = 13.949974, from.ss, num.days = as.numeric(difftime(time1 = to.ss, time2 = from.ss, units = "days", tz = "UTC")), timezone = "UTC")
nightpols <- getNighttimePols(sunrise.sunset.t)


