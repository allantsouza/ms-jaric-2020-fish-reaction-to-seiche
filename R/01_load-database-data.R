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

write_csv(x = hobo_data, file = here("data", "raw", "db", "hobo_data.csv"))





# Interpolated temperatures -----------------------------------------------

sql_query_hobodeployements <- sql_query_temperatures <- paste0(
  "SELECT 
  lp_pos_name as location,
  hl_logger_sn logger_sn,
  hde_depth deployment_depth,
  hde_retrieval_utc retrieval_ts,
  hde_deployment_utc deployment_ts
  FROM at_macfish.hobodeployments
  WHERE lp_pos_name IN ('East', 'West') AND
  (hde_deployment_utc BETWEEN  '", DATE_RANGE[1], "' AND '", DATE_RANGE[2], "' OR hde_retrieval_utc BETWEEN  '", DATE_RANGE[1], "' AND '", DATE_RANGE[2], "')
  ;")

hobodeployements <- data.table(dbGetQuery(con, sql_query_hobodeployements, stringsAsFactors = F))

write_csv(x = hobodeployements, file = here("data", "raw", "db", "hobodeployements.csv"))


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
  st_write(here("data/raw/db/logger_positions.shp"), append = F)


# Lake shape polygon
st_read(con, query = "SELECT lt_shoreline_pol_smooth FROM at_macfish.laketable WHERE lt_lake = 'Chabarovice'") %>%
  st_write(dsn = here("data/raw/db/lake_shape.shp"), append = F)

# Lake depth raster
depth_raster <- rpostgis::pgGetRast(con, c("at_macfish", "depthrast"), clauses = "WHERE lt_lake = 'Chabarovice'")
raster::writeRaster(x = depth_raster, filename = here("data/raw/db/lake_raster.tif"), format = "GTiff", overwrite = T)



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

write_csv(x = wind_data, file = here("data", "raw", "db", "wind_data.csv"))



# Detection data ----------------------------------------------------------
fish_raw <- read_csv(file = here("data/raw/fishIDs.csv"), col_types = "ccdi")

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
  if(nrow(positions) > 0){
    write_csv(x = positions, file = here("data", "raw", "db", "fish", "positions", paste0(tag_sns[i], ".csv")))
  }
}

for(i in 1:length(tag_sns)){
  sql_query_detections <-paste0("
  SELECT dd_id as det_id,
         dd_tag_sn as tag_sn,
         dd_depth as det_depth,
         dd_timestamp_utc as det_ts
  FROM at_macfish.detsdepth
  WHERE fishvalid AND
        dd_tag_sn ='", tag_sns[i], "' AND
        dd_timestamp_utc BETWEEN '", DATE_RANGE[1], "' AND '", DATE_RANGE[2], "'
  ")
  
  detections <-dbGetQuery(con, sql_query_detections, stringsAsFactors = F)
  if(nrow(detections) > 0){
    write_csv(x = detections, file = here("data", "raw", "db", "fish", "detections", paste0(tag_sns[i], ".csv")))
  }
}
