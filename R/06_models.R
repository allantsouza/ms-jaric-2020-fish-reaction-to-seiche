# Loading the data ####
fish_raw <- read_csv(file = "data/raw/fishIDs.csv", col_types = "ccdc") %>%
  mutate(data_path = here("data/products/fish",paste0(tag_sn, ".csv")))

detections <- fish_raw %>%
  filter(file.exists(data_path)) %>%
#  slice(1) %>%   # slicing takes only one factor per level ID
  pull(data_path) %>%
  # read in all the files, appending the path before the filename
  map(~ read_csv(.)) %>% 
  reduce(rbind) %>%
  inner_join(fish_raw[,c("tag_sn","fishid", "species")]) %>%
  dplyr::select(
         fishid,
         species,
         is_valid_seiche,
         diel_period,
         det_depth,
         dets_ts,
         amplitude = det_therm_deviation_crit,
         seasonal_depth = det_location_therm_depth_smoothed_crit,
         mean_gradient = det_therm_gradient_crit
         )

#Transfom variables
detections$dets_ts<-as.integer(detections$dets_ts)  # need to convert time to integer
detections <- detections %>% 
  mutate_at(c("mean_gradient",
              "seasonal_depth",
              "amplitude",
 #             "dets_ts" # time doesn't need to be scaled
  ),
  .funs = ~ scale(.)) #scale predictors only

# All models use only one formula
global_model_formula <- formula(
  det_depth ~
    s(seasonal_depth, k = 10, bs = 'cr') +
    s(amplitude, k = 10, bs = 'cr')+
    s(mean_gradient, k = 10, bs = 'cr') +   # add missing line and remove repeated random smooth for mean_gradient
    s(seasonal_depth, fishid,  bs="fs", m=1) +
    s(amplitude, fishid,  bs="fs", m=1) +
    s(mean_gradient, fishid,  bs="fs", m=1) +
    s(dets_ts, k = 10, bs = 'cr') +         # add smoothing parameter to time
    s(fishid, dets_ts, bs="fs", m=1)
)


#GAMMs with autocorrelation structure ####
##Pike - day ####

#Set pike-day data (from here on to simplify model formula; otherwise it's a mess)
data_pike_day = detections %>%
                             filter(species == "pike" &
                                      diel_period == 'day' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))

#Running the model without autocorrelation to estimate the rho value for the model with autocorrelation
tic('Model run')
mdl_pike_day_simple <- bam(formula = det_depth ~
                             s(lake_therm_thickness_smoothed, k = 100, bs = 'cr') +
                             s(det_therm_strength, k = 100, bs = 'cr') +
                             s(lake_therm_depth_smoothed_center, k = 100, bs = 'cr') +
                             s(det_therm_deviation_center, k = 100, bs = 'cr')+
                             s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                             s(det_therm_strength, fishid,  bs="fs", m=1) +
                             s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) + 
                             s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                             s(dets_ts, k = 100, bs = 'cr') +
                             s(fishid, dets_ts, bs="fs", m=1),
                           data = data_pike_day,  # replace by dataset name
                           family = 'gaussian',
                           nthreads = 10, cluster = 10, gc.level = 0)
toc()

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_pike_day_simple, plot = TRUE)

#Model with autocorrelation
tic('Model run takes')
mld_gamm_pike_day <- bam(formula = det_depth ~
                           s(lake_therm_thickness_smoothed, k = 100, bs = 'cr') +
                           s(det_therm_strength, k = 100, bs = 'cr') +
                           s(lake_therm_depth_smoothed_center, k = 100, bs = 'cr') +
                           s(det_therm_deviation_center, k = 100, bs = 'cr')+
                           s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                           s(det_therm_strength, fishid,  bs="fs", m=1) +
                           s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                           s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                           s(dets_ts, k = 100, bs = 'cr') +
                           s(fishid, dets_ts, bs="fs", m=1),
                         data = data_pike_day,
                         family = 'gaussian',
                         nthreads = 10, cluster = 10, gc.level = 0,
                         AR.start = startindex, rho = rho_start_value) #autocorrelation

toc()

#Model summary
summary.gam(mld_gamm_pike_day)               
#graphical visualization of the model
plot(mld_gamm_pike_day, shade = TRUE, pages = 1, scale = 0)
#Checking the autocorrelation (split by fishid)
acf_resid(model = mld_gamm_pike_day, split_pred = 'fishid')
par(mfrow = c(2,2))
#Checking the model residuals
itsadug::check_resid(model = mld_gamm_pike_day, split_pred = 'fishid')
dev.off()
#Running the model diagnostics
itsadug::diagnostics(model = mld_gamm_pike_day) #it takes a while
#Report stats
itsadug::report_stats(mld_gamm_pike_day)
tab_model(mld_gamm_pike_day, show.ci = F)

##Pike - night ####
#Running the model without autocorrelation to estimate the rho value for the model with autocorrelation
tic('Model run')
mdl_pike_night_simple <- bam(formula = det_depth ~
                             s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                             s(det_therm_strength, k = 10, bs = 'cr') +
                             s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                             s(det_therm_deviation_center, k = 10, bs = 'cr')+
                             s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                             s(det_therm_strength, fishid,  bs="fs", m=1) +
                             s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                             s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                             s(dets_ts) +
                             s(fishid, dets_ts, bs="fs", m=1),
                           data = detections %>%
                             filter(species == "pike" &
                                      diel_period == 'night' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid)),
                           family = 'gaussian',
                           nthreads = 10, cluster = 10, gc.level = 0)
toc()

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_pike_night_simple, plot = TRUE)

#Model with autocorrelation
tic('Model run takes')
mld_gamm_pike_night <- bam(formula = det_depth ~
                             s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                             s(det_therm_strength, k = 10, bs = 'cr') +
                             s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                             s(det_therm_deviation_center, k = 10, bs = 'cr')+
                             s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                             s(det_therm_strength, fishid,  bs="fs", m=1) +
                             s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                             s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                             s(dets_ts) +
                             s(fishid, dets_ts, bs="fs", m=1),
                           data = detections %>%
                             filter(species == "pike" &
                                      diel_period == 'day' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid)),
                           family = 'gaussian',
                           nthreads = 10, cluster = 10, gc.level = 0,
                           AR.start = startindex, rho = rho_start_value) #autocorrelation
toc()

#Model summary
summary.gam(mld_gamm_pike_night)
#graphical visualization of the model
plot(mld_gamm_pike_night, shade = TRUE, pages = 1, scale = 0)
#Checking the autocorrelation (split by fishid)
acf_resid(model = mld_gamm_pike_night, split_pred = 'fishid')
par(mfrow = c(2,2))
#Checking the model residuals
itsadug::check_resid(model = mld_gamm_pike_night, split_pred = 'fishid')
dev.off()
#Running the model diagnostics
itsadug::diagnostics(model = mld_gamm_pike_night) #it takes a while
#Report stats
itsadug::report_stats(mld_gamm_pike_night)
tab_model(mld_gamm_pike_night, show.ci = F)

##Wels - day ####
#Running the model without autocorrelation to estimate the rho value for the model with autocorrelation
tic('Model run')
mdl_wels_day_simple <- bam(formula = det_depth ~
                             s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                             s(det_therm_strength, k = 10, bs = 'cr') +
                             s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                             s(det_therm_deviation_center, k = 10, bs = 'cr')+
                             s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                             s(det_therm_strength, fishid,  bs="fs", m=1) +
                             s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                             s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                             s(dets_ts) +
                             s(fishid, dets_ts, bs="fs", m=1),
                           data = detections %>%
                             filter(species == "wels" &
                                      diel_period == 'day' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid)),
                           family = 'gaussian',
                           nthreads = 10, cluster = 10, gc.level = 0)
toc()

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_wels_day_simple, plot = TRUE)

#Model with autocorrelation
tic('Model run takes')
mld_gamm_wels_day <- bam(formula = det_depth ~
                           s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                           s(det_therm_strength, k = 10, bs = 'cr') +
                           s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                           s(det_therm_deviation_center, k = 10, bs = 'cr')+
                           s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                           s(det_therm_strength, fishid,  bs="fs", m=1) +
                           s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                           s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                           s(dets_ts) +
                           s(fishid, dets_ts, bs="fs", m=1),
                         data = detections %>%
                           filter(species == "wels" &
                                    diel_period == 'day' &
                                    is_valid_seiche == TRUE) %>%
                           mutate(startindex = ifelse(
                             test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                           mutate(fishid = as_factor(fishid)),
                         family = 'gaussian',
                         nthreads = 10, cluster = 10, gc.level = 0,
                         AR.start = startindex, rho = rho_start_value) #autocorrelation

toc()

#Model summary
summary.gam(mld_gamm_wels_day)
#graphical visualization of the model
plot(mld_gamm_wels_day, shade = TRUE, pages = 1, scale = 0)
#Checking the autocorrelation (split by fishid)
acf_resid(model = mld_gamm_wels_day, split_pred = 'fishid')
par(mfrow = c(2,2))
#Checking the model residuals
itsadug::check_resid(model = mld_gamm_wels_day, split_pred = 'fishid')
dev.off()
#Running the model diagnostics
itsadug::diagnostics(model = mld_gamm_wels_day) #it takes a while
#Report stats
itsadug::report_stats(mld_gamm_wels_day)
tab_model(mld_gamm_wels_day, show.ci = F)

##Wels - night ####
#Running the model without autocorrelation to estimate the rho value for the model with autocorrelation
tic('Model run')
mdl_wels_night_simple <- bam(formula = det_depth ~
                               s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                               s(det_therm_strength, k = 10, bs = 'cr') +
                               s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                               s(det_therm_deviation_center, k = 10, bs = 'cr')+
                               s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                               s(det_therm_strength, fishid,  bs="fs", m=1) +
                               s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                               s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                               s(dets_ts) +
                               s(fishid, dets_ts, bs="fs", m=1),
                             data = detections %>%
                               filter(species == "wels" &
                                        diel_period == 'night' &
                                        is_valid_seiche == TRUE) %>%
                               mutate(startindex = ifelse(
                                 test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                               mutate(fishid = as_factor(fishid)),
                             family = 'gaussian',
                             nthreads = 10, cluster = 10, gc.level = 0)
toc()

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_wels_night_simple, plot = TRUE)

#Model with autocorrelation
tic('Model run takes')
mld_gamm_wels_night <- bam(formula = det_depth ~
                             s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                             s(det_therm_strength, k = 10, bs = 'cr') +
                             s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                             s(det_therm_deviation_center, k = 10, bs = 'cr')+
                             s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                             s(det_therm_strength, fishid,  bs="fs", m=1) +
                             s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                             s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                             s(dets_ts) +
                             s(fishid, dets_ts, bs="fs", m=1),
                           data = detections %>%
                             filter(species == "wels" &
                                      diel_period == 'day' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid)),
                           family = 'gaussian',
                           nthreads = 10, cluster = 10, gc.level = 0,
                           AR.start = startindex, rho = rho_start_value) #autocorrelation
toc()

#Model summary
summary.gam(mld_gamm_wels_night)
#graphical visualization of the model
plot(mld_gamm_wels_night, shade = TRUE, pages = 1, scale = 0)
#Checking the autocorrelation (split by fishid)
acf_resid(model = mld_gamm_wels_night, split_pred = 'fishid')
par(mfrow = c(2,2))
#Checking the model residuals
itsadug::check_resid(model = mld_gamm_wels_night, split_pred = 'fishid')
dev.off()
#Running the model diagnostics
itsadug::diagnostics(model = mld_gamm_wels_night) #it takes a while
#Report stats
itsadug::report_stats(mld_gamm_wels_night)
tab_model(mld_gamm_wels_night, show.ci = F)

##Tench - day ####
#Running the model without autocorrelation to estimate the rho value for the model with autocorrelation
tic('Model run')
mdl_tench_day_simple <- bam(formula = det_depth ~
                             s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                             s(det_therm_strength, k = 10, bs = 'cr') +
                             s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                             s(det_therm_deviation_center, k = 10, bs = 'cr')+
                             s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                             s(det_therm_strength, fishid,  bs="fs", m=1) +
                             s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                             s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                             s(dets_ts) +
                             s(fishid, dets_ts, bs="fs", m=1),
                           data = detections %>%
                             filter(species == "tench" &
                                      diel_period == 'day' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid)),
                           family = 'gaussian',
                           nthreads = 10, cluster = 10, gc.level = 0)
toc()

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_tench_day_simple, plot = TRUE)

#Model with autocorrelation
tic('Model run takes')
mld_gamm_tench_day <- bam(formula = det_depth ~
                           s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                           s(det_therm_strength, k = 10, bs = 'cr') +
                           s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                           s(det_therm_deviation_center, k = 10, bs = 'cr')+
                           s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                           s(det_therm_strength, fishid,  bs="fs", m=1) +
                           s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                           s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                           s(dets_ts) +
                           s(fishid, dets_ts, bs="fs", m=1),
                         data = detections %>%
                           filter(species == "tench" &
                                    diel_period == 'day' &
                                    is_valid_seiche == TRUE) %>%
                           mutate(startindex = ifelse(
                             test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                           mutate(fishid = as_factor(fishid)),
                         family = 'gaussian',
                         nthreads = 10, cluster = 10, gc.level = 0,
                         AR.start = startindex, rho = rho_start_value) #autocorrelation

toc()

#Model summary
summary.gam(mld_gamm_tench_day)
#graphical visualization of the model
plot(mld_gamm_tench_day, shade = TRUE, pages = 1, scale = 0)
#Checking the autocorrelation (split by fishid)
acf_resid(model = mld_gamm_tench_day, split_pred = 'fishid')
par(mfrow = c(2,2))
#Checking the model residuals
itsadug::check_resid(model = mld_gamm_tench_day, split_pred = 'fishid')
dev.off()
#Running the model diagnostics
itsadug::diagnostics(model = mld_gamm_tench_day) #it takes a while
#Report stats
itsadug::report_stats(mld_gamm_tench_day)
tab_model(mld_gamm_tench_day, show.ci = F)

##Tench - night ####
#Running the model without autocorrelation to estimate the rho value for the model with autocorrelation
tic('Model run')
mdl_tench_night_simple <- bam(formula = det_depth ~
                               s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                               s(det_therm_strength, k = 10, bs = 'cr') +
                               s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                               s(det_therm_deviation_center, k = 10, bs = 'cr')+
                               s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                               s(det_therm_strength, fishid,  bs="fs", m=1) +
                               s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                               s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                               s(dets_ts) +
                               s(fishid, dets_ts, bs="fs", m=1),
                             data = detections %>%
                               filter(species == "tench" &
                                        diel_period == 'night' &
                                        is_valid_seiche == TRUE) %>%
                               mutate(startindex = ifelse(
                                 test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                               mutate(fishid = as_factor(fishid)),
                             family = 'gaussian',
                             nthreads = 10, cluster = 10, gc.level = 0)
toc()

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_tench_night_simple, plot = TRUE)

#Model with autocorrelation
tic('Model run takes')
mld_gamm_tench_night <- bam(formula = det_depth ~
                             s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                             s(det_therm_strength, k = 10, bs = 'cr') +
                             s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                             s(det_therm_deviation_center, k = 10, bs = 'cr')+
                             s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                             s(det_therm_strength, fishid,  bs="fs", m=1) +
                             s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                             s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                             s(dets_ts) +
                             s(fishid, dets_ts, bs="fs", m=1),
                           data = detections %>%
                             filter(species == "tench" &
                                      diel_period == 'day' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid)),
                           family = 'gaussian',
                           nthreads = 10, cluster = 10, gc.level = 0,
                           AR.start = startindex, rho = rho_start_value) #autocorrelation
toc()

#Model summary
summary.gam(mld_gamm_tench_night)
#graphical visualization of the model
plot(mld_gamm_tench_night, shade = TRUE, pages = 1, scale = 0)
#Checking the autocorrelation (split by fishid)
acf_resid(model = mld_gamm_tench_night, split_pred = 'fishid')
par(mfrow = c(2,2))
#Checking the model residuals
itsadug::check_resid(model = mld_gamm_tench_night, split_pred = 'fishid')
dev.off()
#Running the model diagnostics
itsadug::diagnostics(model = mld_gamm_tench_night) #it takes a while
#Report stats
itsadug::report_stats(mld_gamm_tench_night)
tab_model(mld_gamm_tench_night, show.ci = F)

##Rudd - day ####
#Running the model without autocorrelation to estimate the rho value for the model with autocorrelation
tic('Model run')
mdl_rudd_day_simple <- bam(formula = det_depth ~
                              s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                              s(det_therm_strength, k = 10, bs = 'cr') +
                              s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                              s(det_therm_deviation_center, k = 10, bs = 'cr')+
                              s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                              s(det_therm_strength, fishid,  bs="fs", m=1) +
                              s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                              s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                              s(dets_ts) +
                              s(fishid, dets_ts, bs="fs", m=1),
                            data = detections %>%
                              filter(species == "rudd" &
                                       diel_period == 'day' &
                                       is_valid_seiche == TRUE) %>%
                              mutate(startindex = ifelse(
                                test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                              mutate(fishid = as_factor(fishid)),
                            family = 'gaussian',
                            nthreads = 10, cluster = 10, gc.level = 0)
toc()

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_rudd_day_simple, plot = TRUE)

#Model with autocorrelation
tic('Model run takes')
mld_gamm_rudd_day <- bam(formula = det_depth ~
                            s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                            s(det_therm_strength, k = 10, bs = 'cr') +
                            s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                            s(det_therm_deviation_center, k = 10, bs = 'cr')+
                            s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                            s(det_therm_strength, fishid,  bs="fs", m=1) +
                            s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                            s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                            s(dets_ts) +
                            s(fishid, dets_ts, bs="fs", m=1),
                          data = detections %>%
                            filter(species == "rudd" &
                                     diel_period == 'day' &
                                     is_valid_seiche == TRUE) %>%
                            mutate(startindex = ifelse(
                              test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                            mutate(fishid = as_factor(fishid)),
                          family = 'gaussian',
                          nthreads = 10, cluster = 10, gc.level = 0,
                          AR.start = startindex, rho = rho_start_value) #autocorrelation

toc()

#Model summary
summary.gam(mld_gamm_rudd_day)
#graphical visualization of the model
plot(mld_gamm_rudd_day, shade = TRUE, pages = 1, scale = 0)
#Checking the autocorrelation (split by fishid)
acf_resid(model = mld_gamm_rudd_day, split_pred = 'fishid')
par(mfrow = c(2,2))
#Checking the model residuals
itsadug::check_resid(model = mld_gamm_rudd_day, split_pred = 'fishid')
dev.off()
#Running the model diagnostics
itsadug::diagnostics(model = mld_gamm_rudd_day) #it takes a while
#Report stats
itsadug::report_stats(mld_gamm_rudd_day)
tab_model(mld_gamm_rudd_day, show.ci = F)

##Rudd - night ####
#Running the model without autocorrelation to estimate the rho value for the model with autocorrelation
tic('Model run')
mdl_rudd_night_simple <- bam(formula = det_depth ~
                                s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                                s(det_therm_strength, k = 10, bs = 'cr') +
                                s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                                s(det_therm_deviation_center, k = 10, bs = 'cr')+
                                s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                                s(det_therm_strength, fishid,  bs="fs", m=1) +
                                s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                                s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                                s(dets_ts) +
                                s(fishid, dets_ts, bs="fs", m=1),
                              data = detections %>%
                                filter(species == "rudd" &
                                         diel_period == 'night' &
                                         is_valid_seiche == TRUE) %>%
                                mutate(startindex = ifelse(
                                  test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                                mutate(fishid = as_factor(fishid)),
                              family = 'gaussian',
                              nthreads = 10, cluster = 10, gc.level = 0)
toc()

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_rudd_night_simple, plot = TRUE)

#Model with autocorrelation
tic('Model run takes')
mld_gamm_rudd_night <- bam(formula = det_depth ~
                              s(lake_therm_thickness_smoothed, k = 10, bs = 'cr') +
                              s(det_therm_strength, k = 10, bs = 'cr') +
                              s(lake_therm_depth_smoothed_center, k = 10, bs = 'cr') +
                              s(det_therm_deviation_center, k = 10, bs = 'cr')+
                              s(lake_therm_thickness_smoothed, fishid,  bs="fs", m=1) +
                              s(det_therm_strength, fishid,  bs="fs", m=1) +
                              s(lake_therm_depth_smoothed_center, fishid,  bs="fs", m=1) +
                              s(det_therm_deviation_center, fishid,  bs="fs", m=1) +
                              s(dets_ts) +
                              s(fishid, dets_ts, bs="fs", m=1),
                            data = detections %>%
                              filter(species == "rudd" &
                                       diel_period == 'day' &
                                       is_valid_seiche == TRUE) %>%
                              mutate(startindex = ifelse(
                                test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                              mutate(fishid = as_factor(fishid)),
                            family = 'gaussian',
                            nthreads = 10, cluster = 10, gc.level = 0,
                            AR.start = startindex, rho = rho_start_value) #autocorrelation
toc()

#Model summary
summary.gam(mld_gamm_rudd_night)
#graphical visualization of the model
plot(mld_gamm_rudd_night, shade = TRUE, pages = 1, scale = 0)
#Checking the autocorrelation (split by fishid)
acf_resid(model = mld_gamm_rudd_night, split_pred = 'fishid')
par(mfrow = c(2,2))
#Checking the model residuals
itsadug::check_resid(model = mld_gamm_rudd_night, split_pred = 'fishid')
dev.off()
#Running the model diagnostics
itsadug::diagnostics(model = mld_gamm_rudd_night) #it takes a while
#Report stats
itsadug::report_stats(mld_gamm_rudd_night)
tab_model(mld_gamm_rudd_night, show.ci = F)
