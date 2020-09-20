# Loading all environmental data from database


# Hobo data ---------------------------------------------------------------

sql_query_hobo_data <- paste0(
  "SELECT
    a.hl_logger_sn loggerSN,
    hd_timestamp_utc ts,
    hd_temperature temperature,
    hd_light_lux light,
    hde_depth depth,
    lp_pos_name as location,
    lt_lake lake,
    hde_deploymentid deployment_id,
    hde_retrieval_utc ts_retrieval,
    hde_deployment_utc  as ts_deployment
  FROM (
      SELECT lt_lake,
        lp_pos_order,
        hde_deploymentid,
        hl_logger_sn,
        hde_depth,
        x.lp_pos_name,
        hde_retrieval_utc,
        hde_deployment_utc
      FROM at_macfish.hobodeployments x 
      INNER JOIN at_macfish.loggerpos y 
      ON x.lp_pos_name = y.lp_pos_name
      WHERE x.lp_pos_name IN ('East', 'West')
      ) a 
  INNER JOIN (SELECT *
      FROM at_macfish.hobodata
      WHERE hd_timestamp_utc BETWEEN '", DATE_RANGE[1], "' AND '", DATE_RANGE[2], "' AND
      hd_valid
      ) l 
  ON l.hl_logger_sn = a.hl_logger_sn AND
  hd_timestamp_utc BETWEEN hde_deployment_utc AND hde_retrieval_utc"
)

hobo_data <- dbGetQuery(conn = con, sql_query_hobo_data)
tz(hobo_data$ts) <- "UTC"
hobo_data <- data.table(hobo_data)

#filter out logs with NA or inf temperature
hobo_data <- hobo_data[is.finite(temperature)]
#round timestamp to have 5 minutes interval - each logger recorded 1 log per 5 minute
#add one second since one data logger was logging exactly in 5 min intervals
hobo_data$interval <- lubridate::round_date(hobo_data$ts, "5 minutes")

#for each log, define its order for each 5 min interval]
setkey(hobo_data, location, depth)
hobo_data[, order := 1:.N, by =.(location, as.numeric(interval))]

write_csv(x = hobo_data, path = here("data", "raw", "db", "hobo_data.csv"))





# Interpolated temperatures -----------------------------------------------

# load interpolated temperatures from the desired range of Dates
# TODO: interpolate hobo_data locally so it is clear how that was done
sql_query_temperatures <- paste0(
  "SELECT 
  hte_pos_name as location,
  hte_timestamp_utc ts,
  hte_depth depth,
  hte_temperature temperature  
  FROM at_macfish.hobo_temperatures
  WHERE hte_pos_name IN ('East', 'West') AND
  hte_timestamp_utc BETWEEN  '", DATE_RANGE[1], "' AND '", DATE_RANGE[2], "';") 

temperatures <- data.table(dbGetQuery(con, sql_query_temperatures, stringsAsFactors = F))
tz(temperatures$ts) <- "UTC"
#Since there are some time frames in which the temperature increases over depth for some steps ( location == "West" & ts == "2015-06-23 02:35:00")
setkey(temperatures, location, ts, depth)
temperatures[, non_monotonic_decrease := any(diff(temperature) > 0), by = .(location, ts)]
temperatures[non_monotonic_decrease == T] # TODO:
#remove those steps and interpolate linearly created gaps
temperatures[, cum.min.temp := cummin(temperature), by = .(location, ts)]  
temperatures[cum.min.temp < temperature, temperature := cum.min.temp]
#since you are interested in last occurence of that temperature, you can exclude those above
temperatures_monotonic <- temperatures[, .(depth = max(depth)) , by = .(location, ts, temperature)]

write_csv(x = temperatures_monotonic , path = here("data", "raw", "db", "temperature_data.csv"))



# Spatial objects ---------------------------------------------------------

sql_query_loggerpos <- "
SELECT 
  lp_pos_name as location,
  st_x(lp_geom) easting,
  st_y(lp_geom) northing,
  lt_lake lake,
  lp_geom,
  lp_pos_order location_order
FROM at_macfish.loggerpos 
WHERE lt_lake = 'Chabarovice'
"

st_read(con, query = sql_query_loggerpos) %>% 
  st_write("data/raw/db/logger_positions.shp")


# Lake shape polygon
st_read(con, query = "SELECT lt_shoreline_pol_smooth FROM at_macfish.laketable WHERE lt_lake = 'Chabarovice'") %>%
  st_write(obj = lake_shp, dsn = "data/raw/db/lake_shape.shp" )



# Wind data ---------------------------------------------------------------

#cumulative wind run over 30 minutes
windrun <- data.table(dbGetQuery(conn = con, statement = "SELECT * FROM at_macfish.meteodata_view where mds_dataset_id IN (22)"))
#mean direction of wind in degrees of azimuth
winddirection <- data.table(dbGetQuery(conn = con, statement = "SELECT * FROM at_macfish.meteodata_view where mds_dataset_id IN (32)"))
#mean windspeed
windspeed <- data.table(dbGetQuery(conn = con, statement = "SELECT * FROM at_macfish.meteodata_view where mds_dataset_id IN (20)"))

wind_data <- rbind(windrun, windspeed, winddirection)
wind_data <- rename(wind_data, 
       lake = ms_lake, 
       ts = md_timestamp_utc,
       dataset_id = mds_dataset_id)

write_csv(x = wind_data, path = here("data", "raw", "db", "wind_data.csv"))



# Detection data ----------------------------------------------------------
fish_raw <- read_csv(file = "data/raw/fishIDs.csv", col_types = "ccd")

tag_sns <- fish_raw$tag_sn

for(i in 1:length(tag_sns)){
  sql_query_positions <-paste0("
  SELECT up_tag_sn, up_id, up_disttoshore, up_timestamp_utc, up_depth, up_lake lake, ST_X(up_geom) x, ST_Y(up_geom) y 
     FROM at_macfish.umap_pos_powfilter 
     WHERE up_validpos AND
           fishvalid AND
           up_tag_sn ='", tag_sns[i], "' AND
           up_timestamp_utc BETWEEN '", DATE_RANGE[1], "' AND '", DATE_RANGE[2], "'
           ")
  positions <- dbGetQuery(con, sql_query_positions, stringsAsFactors = F)
  write_csv(x = positions, path = here("data", "raw", "db", "fish", "positions", paste0(tag_sns[i], ".csv")))
  
  
  sql_query_detections <-paste0("
  SELECT dd_tag_sn as tag_sn,
         dd_depth as det_depth,
         dd_timestamp_utc as det_ts
  FROM at_macfish.detsdepth
  WHERE fishvalid AND
        dd_tag_sn ='", tag_sns[i], "' AND
        dd_timestamp_utc BETWEEN '", DATE_RANGE[1], "' AND '", DATE_RANGE[2], "'
  ")
  detections <-dbGetQuery(con, sql_query_detections, stringsAsFactors = F)
  
  write_csv(x = detections, path = here("data", "raw", "db", "fish", "detections", paste0(tag_sns[i], ".csv")))
  
}
