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

# All models use this formula as a starting point: Models are re-fitted with increasing number of knots to check changes in effective degrees of freedom (edf).
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


################################################################################ Pike - day #####################################################################################

# Set pike-day data (from here on to simplify model formula; otherwise it's a mess)
data_pike_day = detections %>%
                             filter(species == "pike" &
                                      diel_period == 'day' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))


#######################################
### CHOOSING THE RIGHT DISTRIBUTION ###
#######################################
# Some dataset seems right or left skewed, others appear to require mixed distributions (e.g. Tweedie).
# Fitting models with Gamma or inverse Gaussian is very time-consuming, more so if we need to adjust the basis dimensions and re-fit models.
# I introduce the library "fitdistrplus" to help decide which distribution fits better to each of the subsets analyzed.

library(fitdistrplus)
library(logspline)

det_depth <- data_pike_day$det_depth

# Kurtosis vs. squared skewness 
par(mfrow = c(2,2))
descdist(det_depth, discrete = FALSE)   # normal, gamma and weibull are possible distributions  

# Comparison of various distribution fits 
fitW <- fitdist(det_depth, "weibull")
fitg <- fitdist(det_depth, "gamma")
fitln <- fitdist(det_depth, "lnorm")
fitn <- fitdist(det_depth, "norm")
fitexp <- fitdist(det_depth, "exp")
fitunif <- fitdist(det_depth, "unif")
summary(fitW)
summary(fitg)
summary(fitln)
summary(fitn)
summary(fitexp)
summary(fitunif)
cdfcomp(list(fitW, fitg, fitln, fitn, fitexp, fitunif), legendtext=c("Weibull", "gamma", "lognormal", "normal", "exponential", "uniform"))
denscomp(list(fitW, fitg, fitln, fitn, fitexp, fitunif), legendtext=c("Weibull", "gamma", "lognormal", "normal", "exponential", "uniform"))
qqcomp(list(fitW, fitg, fitln, fitn, fitexp, fitunif), legendtext=c("Weibull", "gamma", "lognormal", "normal", "exponential", "uniform"))
ppcomp(list(fitW, fitg, fitln, fitn, fitexp, fitunif), legendtext=c("Weibull", "gamma", "lognormal", "normal", "exponential", "uniform"))
gofstat(list(fitW, fitg, fitln, fitn, fitexp, fitunif), fitnames=c("Weibull", "gamma", "lognormal", "normal", "exponential", "uniform"))

Goodness-of-fit statistics
                                  Weibull        gamma    lognormal       normal  exponential      uniform
Kolmogorov-Smirnov statistic 1.804174e-01 2.228361e-01 2.430853e-01    0.1616309 3.001863e-01 4.980728e-01
Cramer-von Mises statistic   1.700514e+03 2.671305e+03 3.436801e+03 1272.9246927 7.458845e+03 1.570884e+04
Anderson-Darling statistic   1.011722e+04 1.413371e+04 1.820985e+04 7219.1229215 3.758546e+04          Inf

Goodness-of-fit criteria
                               Weibull   gamma lognormal  normal exponential uniform
Akaikes Information Criterion  1232376 1290480   1357686 1218099     1452296 1500518    
Bayesian Information Criterion 1232396 1290501   1357706 1218120     1452306 1500539      # In this case, normal seems a reasonable option


#########################################################
### MODEL SELECTION BY ADJUSTMENT OF BASIS DIMENSIONS ###
#########################################################
# The same base model is fitted with different number of knots (k) (10, 100, 200, 500 and 1000[1440 per day]) to test wether there are changes in the edf of the predictors.
# When basis dimensions are estabilised between models we select the one where changes start to be insignificant. 
# Next, extract the edf for each smooth in the selected model and round them.
# Re-fit the model with fixed edf (without penalization).
# Check if the smooths of some predictors have a linear relationship and in that case, drop the smoothing function and re-fit the model.
# Compare the full-smoothing and the simplified models.
# To summarize, only the model formulas of the selected model (before and after the edf fix) are shown.
# This procedure has proven to be optimal for reducing overdispersion (for example, from overfitting), increasing both R^2 and the deviance explained and to get more accurate p-values. 
# In fact, some predictor terms changed completely their significance before and after the adjustments or they had unreliable low p-values (other than the random smooths). 
# In this case (pike-day), k=100 is our choice. The model is subsequently re-fitted by fixing edf and dropping one smooth function (mean_gradient). 

# 1.1. Running the model without autocorrelation to estimate the rho value for the model with autocorrelation
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

# 1.2. Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_pike_day_simple, plot = TRUE)

# 1.3. Model with autocorrelation
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

# Checking the autocorrelation (split by fishid)
acf_resid(model = mld_gamm_pike_day, split_pred = 'fishid')
par(mfrow = c(2,2))

# Model summary
summary.gam(mld_gamm_pike_day)

# Extract R^2 
summary(mld_gamm_pike_day)$r.sq   
[1] 0.745     # R^2 is not bad but can improve

# Check overdispersion (phi)
sum(resid(mld_gamm_pike_day, type = "pearson")^2) / df.residual(mld_gamm_pike_day)
[1] 2.380257  # too high (phi > 1.5)

# 2.1. Extract and round edf of the smooth terms (estimated with penalization) to fix them
f.df <- round(summary(mld_gamm_pike_day_cr_k100)$edf)+1  # get edf per smooth
f.df <- pmax(f.df,3)                                     # minimum basis dimension is 3

# 2.2. Re-fit the model with fixed edf for each smooth term (unpenalized using fx=TRUE)
mld_gamm_pike_day_edf <- bam(formula =
                         det_depth ~
                         s(seasonal_depth, k=f.df[1], fx=TRUE) +
                         s(amplitude, k=f.df[2], fx=TRUE)+
                         s(mean_gradient, k=f.df[3], fx=TRUE) +
                         s(seasonal_depth, fishid,  k=f.df[4], bs="fs", m=1) +
                         s(amplitude, fishid, k=f.df[5], bs="fs", m=1) +
                         s(mean_gradient, fishid, k=f.df[6], bs="fs", m=1) +
                         s(dets_ts, k=f.df[7], fx=TRUE)+
                         s(fishid, dets_ts, k=f.df[8], bs="fs", m=1),
                         data = data_pike_day,
                         family = 'gaussian',
                         nthreads = 10, cluster = 10, gc.level = 0,
                         AR.start = startindex, rho = rho_start_value)


# 2.3. Compare models 
AIC(mld_gamm_pike_day_cr_edf,mld_gamm_pike_day_cr_k50, mld_gamm_pike_day_cr_k50_edf, mld_gamm_pike_day, mld_gamm_pike_day_edf)

                                df      AIC
mld_gamm_pike_day_cr_edf      2244.0729 214831.3
mld_gamm_pike_day_cr_k50       441.6427 239347.4
mld_gamm_pike_day_cr_k50_edf  2248.3288 213448.5    # The 50k model with fixed edf is very close to 100k
mld_gamm_pike_day              521.2288 237957.6
mld_gamm_pike_day_edf         2281.9432 213072.3    # The 100k model with fixed edf has the lowest AIC

# Is the 100k model with fixed edf preferred to 50k?
compareML(mld_gamm_pike_day_cr_k50_edf, mld_gamm_pike_day_edf)  # 100k fx still has lower fREML score and AIC 
# We can also compare models with and without fixed edf
compareML(mld_gamm_pike_day, mld_gamm_pike_day_edf)             # The fixed model is better

# Model summary
summary.gam(mld_gamm_pike_day_edf)

# Extract R^2 
summary(mld_gamm_pike_day_edf)$r.sq   
[1] 0.88        # R^2 has increased by ~15%/20% with respect to 100-k-non-fx model and default (10k)-non-fx model, respectively

# Check overdispersion (phi)
sum(resid(mld_gamm_pike_day_edf, type = "pearson")^2) / df.residual(mld_gamm_pike_day_edf)
[1] 1.13283    # phi has decreased more than the unit (phi < 1.5)

# Check residuals
qq.gam(mld_gamm_pike_day_edf)             # Much better residuals than previous models
gam.check(mld_gamm_pike_day_edf)          # The k' index seems quite good for all terms

# Plot non-parametric predictors and random smooth effects
par(mfrow = c(2,4))
plot(mld_gamm_pike_day_edf, select = 1, shade = TRUE, scale = 0, seWithMean = TRUE)  # seasonal depth effect
plot(mld_gamm_pike_day_edf, select = 2, shade = TRUE, scale = 0, seWithMean = TRUE)  # amplitude effect
plot(mld_gamm_pike_day_edf, select = 3, shade = TRUE, scale = 0, seWithMean = TRUE)  # mean_gradient effect - Note the linear relationship with det_depth
plot(mld_gamm_pike_day_edf, select = 4, shade = TRUE, scale = 0, seWithMean = TRUE)  # slope of seasonal_depth for each fish
plot(mld_gamm_pike_day_edf, select = 5, shade = TRUE, scale = 0, seWithMean = TRUE)  # slope of amplitude for each fish
plot(mld_gamm_pike_day_edf, select = 6, shade = TRUE, scale = 0, seWithMean = TRUE)  # slope of mean_gradient for each fish
plot(mld_gamm_pike_day_edf, select = 6, shade = TRUE, scale = 0, seWithMean = TRUE)  # time effect
plot(mld_gamm_pike_day_edf, select = 6, shade = TRUE, scale = 0, seWithMean = TRUE)  # slope of time for each fish

# The plots evidence that the smooth for mean_gradient shows linearity.
# Additionally, edf keeps low and close to linear between models independently of the numer of knots, as previously tested so far.
# Hence, we can enter mean_gradient as a linear parametric term.

# 3.1. Re-fit a simplified model equivalent to mld_gamm_pike_day_edf

mld_gamm_pike_day_final <- bam(formula =
                           det_depth ~
                           s(seasonal_depth, k=f.df[1], fx=TRUE) +
                           s(amplitude, k=f.df[2], fx=TRUE)+
                           mean_gradient +
                           s(seasonal_depth, fishid,  k=f.df[4], bs="fs", m=1) +
                           s(amplitude, fishid, k=f.df[5], bs="fs", m=1) +
                           s(mean_gradient, fishid, k=f.df[6], bs="fs", m=1) +
                           s(dets_ts, k=f.df[7], fx=TRUE)+
                           s(fishid, dets_ts, k=f.df[8], bs="fs", m=1),
                           data = data_pike_day,
                           family = 'gaussian',
                           nthreads = 10, cluster = 10, gc.level = 0,
                           AR.start = startindex, rho = rho_start_value)


# 3.2. Check that the full-smoothing and simplified models are identical

AIC(mld_gamm_pike_day_cr_k100_edf, mld_gamm_pike_day_final)
                                    df      AIC
mld_gamm_pike_day_cr_k100_edf 2281.943 213072.3
mld_gamm_pike_day_final       2280.522 213072.6     # There are no differences between models but... 

deviance(mld_gamm_pike_day_edf)
[1] 269419                  
deviance(mld_gamm_pike_day_final)                   # ...deviance still improves a bit after dropping the smooth term for mean_gradient
[1] 269398.3                                       

# Report stats
itsadug::report_stats(mld_gamm_pike_day)
tab_model(mld_gamm_pike_day, show.ci = F)
gamtabs(mld_gamm_pike_day_final, caption="Summaty of mld_gamm_pike_day_final", comment=FALSE, type='html')

################################################################################ Pike - night #####################################################################################

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
