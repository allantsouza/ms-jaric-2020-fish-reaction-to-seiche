# <a name="headindex">Influence of internal seiche dynamics on vertical distribution of fish
## <a name="headindex">Report of model analysis
---

# [Data](#headdata)

# [Final models ](#headafinalmodels)

[![Image name](/outputs/icons/pike_index_1.png)](#headSpeciesPike)
[![Image name](/outputs/icons/wels_index_1.png)](#headSpeciesWels)
[![Image name](/outputs/icons/rudd_index_1.png)](#headSpeciesRudd)
[![Image name](/outputs/icons/tench_index_1.png)](#headSpeciesTench)

# [Final models table ](#headafinaltable)

# [Alternative model structures](#headaltmodels)

# [Issues and further improvements](#headissuesandimprovements)

# [References](#headreferences)

---

# <a name="headdata"></a>Data [:page_facing_up:](#headindex)

Load libraries.

[:books:](https://cran.r-project.org/web/packages/mgcv/index.html)`library(mgcv)`
[:books:](https://cran.r-project.org/web/packages/itsadug/index.html)`library(itsadug)`
[:books:](https://cran.r-project.org/web/packages/gratia/index.html)`library(gratia)`
[:books:](https://cran.r-project.org/web/packages/AICcmodavg/index.html)`library(AICcmodavg)`
[:books:](https://cran.r-project.org/web/packages/knitr/index.html)`library(knitr)`
[:books:](https://cran.r-project.org/web/packages/visreg/index.html)`library(visreg)`
[:books:](https://cran.r-project.org/web/packages/ggplot2/index.html)`library(ggplot2)`
[:books:](https://cran.r-project.org/web/packages/visreg/index.html)`library(visreg)`
[:books:](https://cran.r-project.org/web/packages/ggplot2/index.html)`library(ggplot2)`
[:books:](https://cran.r-project.org/web/packages/fitdistrplus/index.html)`library(fitdistrplus)`
[:books:](https://cran.r-project.org/web/packages/logspline/index.html)`library(logspline)`

Load the whole dataset.
``` r
fish_raw <- read_csv(file = "data/raw/fishIDs.csv", col_types = "ccdc") %>%
  mutate(data_path = here("data/products/fish",paste0(tag_sn, ".csv")))

detections <- fish_raw %>%
  filter(file.exists(data_path)) %>%
  pull(data_path) %>%
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
```

Transform continous covariates (scaling by 2 SD).
``` r

detections$dets_ts<-as.integer(detections$dets_ts)

detections <- detections %>%
  mutate_at(c("mean_gradient",
              "seasonal_depth",
              "amplitude"
  ),
  .funs = ~ scale(.))
```

Summary of the whole dataset.

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:left;">
Description
</th>
<th style="text-align:left;">
Class
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
<span
style="     color: #273746 !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #F2F3F4 !important;">fishid</span>
</td>
<td style="text-align:left;width: 30em; ">
The Fish ID
</td>
<td style="text-align:left;font-style: italic;">
string with 12-19 unique values (depending on data subset)
</td>
</tr>
<tr>
<td style="text-align:left;">
<span
style="     color: #273746 !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #F2F3F4 !important;">diel_period</span>
</td>
<td style="text-align:left;width: 30em; ">
The diel period (day, night)
</td>
<td style="text-align:left;font-style: italic;">
string with 2 unique values
</td>
</tr>
<tr>
<td style="text-align:left;">
<span
style="     color: #273746 !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #F2F3F4 !important;">dets_ts</span>
</td>
<td style="text-align:left;width: 30em; ">
The time series in measurements per minute
</td>
<td style="text-align:left;font-style: italic;">
integer ranging from June 16 to October 15
</td>
</tr>
<tr>
<td style="text-align:left;">
<span
style="     color: #273746 !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #F2F3F4 !important;">seasonal_depth</span>
</td>
<td style="text-align:left;width: 30em; ">
The long-term depth of thermocline with no occurrence of seiche effect
</td>
<td style="text-align:left;font-style: italic;">
numeric ranging from -5 to 5
</td>
</tr>
<tr>
<td style="text-align:left;">
<span
style="     color: #273746 !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #F2F3F4 !important;">amplitude</span>
</td>
<td style="text-align:left;width: 30em; ">
The deviation of actual thermocline depth from seasonal depth
</td>
<td style="text-align:left;font-style: italic;">
numeric ranging from -10 to 10
</td>
</tr>
<tr>
<td style="text-align:left;">
<span
style="     color: #273746 !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #F2F3F4 !important;">mean_gradient</span>
</td>
<td style="text-align:left;width: 30em; ">
The temperature gradient in degree per meter between thermocline top and base
</td>
<td style="text-align:left;font-style: italic;">
numeric ranging from -5 to 12
</td>
</tr>
<tr>
<td style="text-align:left;">
<span
style="     color: #273746 !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #F2F3F4 !important;">det_depth</span>
</td>
<td style="text-align:left;width: 30em; ">
The depth for a fish at a specific time
</td>
<td style="text-align:left;font-style: italic;">
numeric ranging from 0.5 to 30
</td>
</tr>
<tr>
<td style="text-align:left;">
<span
style="     color: #273746 !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: #F2F3F4 !important;">is_valid_seiche</span>
</td>
<td style="text-align:left;width: 30em; ">
Whether the the seiche effect is valid according the lake disbalance of interpolated wide position of thermocline
</td>
<td style="text-align:left;font-style: italic;">
logical
</td>
</tr>
</tbody>
</table>

Generate data subsets for each species and diel periods.

``` r
data_pike_day = detections %>%
                             filter(species == "pike" &
                                      diel_period == 'day' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))

data_pike_night = detections %>%
                             filter(species == "pike" &
                                      diel_period == 'night' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))

data_rudd_day = detections %>%
                             filter(species == "rudd" &
                                      diel_period == 'day' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))

data_rudd_night = detections %>%
                             filter(species == "rudd" &
                                      diel_period == 'night' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))

data_tench_day = detections %>%
                             filter(species == "tench" &
                                      diel_period == 'day' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))

data_tench_night = detections %>%
                             filter(species == "tench" &
                                      diel_period == 'night' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))

data_wels_day = detections %>%
                             filter(species == "wels" &
                                      diel_period == 'day' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))

data_wels_night = detections %>%
                             filter(species == "wels" &
                                      diel_period == 'night' &
                                      is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))
```

Calculate 5-min average fish positions and put them into new datasets
``` r
datasets <- list(data_tench_day,
                 data_tench_night,
          	     data_pike_day,
         	       data_pike_night,
         	       data_wels_day,
         	       data_wels_night,
         	       data_rudd_day,
         	       data_rudd_night)

for(df in datasets){
df_mean_5min <- function(df){
            df <- as.data.frame(df)
            df$dets_ts5min <- floor_date(df$dets_ts, "5 mins")
            df <- df %>% group_by(fishid, dets_ts5min, startindex) %>%
                         summarise_at(vars(det_depth,seasonal_depth,amplitude,mean_gradient),mean)
            df$dets_ts5min<-as.integer(df$dets_ts5min)
            print(df)
}}

data_tench_day_5min <- df_mean_5min(data_tench_day)
data_tench_night_5min <- df_mean_5min(data_tench_night)
data_pike_day_5min <- df_mean_5min(data_pike_day)
data_pike_night_5min <- df_mean_5min(data_pike_night)
data_rudd_day_5min <- df_mean_5min(data_rudd_day)
data_rudd_night_5min <- df_mean_5min(data_rudd_night)
data_wels_day_5min <- df_mean_5min(data_wels_day)
data_wels_night_5min <- df_mean_5min(data_wels_night)
```

# <a name="headafinalmodels"></a> Fit random-intercept-only models [:page_facing_up:](#headindex)

These individual-level models serve to calculate the raw phenotypic variance without controlling for fixed-effects of covariates (model_fishid)
``` r
mdl_tench_day_simple_fishid <- bam(formula = det_depth ~
                                             s(fishid, k = 19, bs = 're'),
                                             data = data_tench_day,
                                             family = 'gaussian', select=TRUE,method="REML",
                                             nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_tench_day_simple_fishid, plot = TRUE)

#Model with autocorrelation
mdl_tench_day_fishid <- bam(formula = det_depth ~
                                      s(fishid, k = 19, bs = 're'),
                                      data = data_tench_day,
                                      family = 'gaussian', select=TRUE,method="REML",
                                      nthreads = 10, cluster = 10, gc.level = 0,
                                      AR.start = startindex, rho = rho_start_value)


mdl_tench_night_simple_fishid <- bam(formula = det_depth ~
                                               s(fishid, k = 19, bs = 're'),
                                               data = data_tench_night_5min,
                                               family = Gamma(link = "log"), select=TRUE,method="REML",
                                               nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_tench_night_simple_fishid, plot = TRUE)

#Model with autocorrelation
mdl_tench_night_fishid <- bam(formula = det_depth ~
                                        s(fishid, k = 19, bs = 're'),
                                        data = data_tench_night_5min,
                                        family = Gamma(link = "log"), select=TRUE,method="REML",
                                        nthreads = 10, cluster = 10, gc.level = 0,
                                        AR.start = startindex, rho = rho_start_value)


mdl_pike_day_simple_fishid <- bam(formula = det_depth ~
                                            s(fishid, k = 19, bs = 're'),
                                            data = data_pike_day_5min,
                                            family = 'gaussian', select=TRUE,method="REML",
                                            nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_pike_day_simple_fishid, plot = TRUE)

#Model with autocorrelation
mdl_pike_day__fishid <- bam(formula = det_depth ~
                                      s(fishid, k = 19, bs = 're'),
                                      data = data_pike_day_5min,
                                      family = 'gaussian', select=TRUE,method="REML",
                                      nthreads = 10, cluster = 10, gc.level = 0,
                                      AR.start = startindex, rho = rho_start_value)



mdl_pike_night_simple_fishid <- bam(formula = det_depth ~
                                              s(fishid, k = 19, bs = 're'),
                                              data = data_pike_night_5min,
                                              family = 'gaussian', select=TRUE,method="REML",
                                              nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_pike_night_simple_fishid, plot = TRUE)

#Model with autocorrelation
mdl_pike_night_fishid <- bam(formula = det_depth ~
                                       s(fishid, k = 19, bs = 're'),
                                       data = data_pike_night_5min,
                                       family = 'gaussian', select=TRUE,method="REML",
                                       nthreads = 10, cluster = 10, gc.level = 0,
                                       AR.start = startindex, rho = rho_start_value)


mdl_wels_day_simple_fishid <- bam(formula = det_depth ~
                                            s(fishid, k = 19, bs = 're'),
                                            data = data_wels_day_5min,
                                            family = Gamma(link = "log"), select=TRUE,method="REML",
                                            nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_wels_day_simple_fishid, plot = TRUE)

#Model with autocorrelation
mdl_wels_day_fishid <- bam(formula = det_depth ~
                                     s(fishid, k = 19, bs = 're'),
                                     data = data_wels_day_5min,
                                     family = Gamma(link = "log"), select=TRUE,method="REML",
                                     nthreads = 10, cluster = 10, gc.level = 0,
                                     AR.start = startindex, rho = rho_start_value)

mdl_wels_night_simple_fishid <- bam(formula = det_depth ~
                                              s(fishid, k = 19, bs = 're'),
                                              data = data_wels_night_5min,
                                              family = Gamma(link = "log"), select=TRUE,method="REML",
                                              nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_wels_night_simple_fishid, plot = TRUE)

#Model with autocorrelation
mdl_wels_night_fishid <- bam(formula = det_depth ~
                                       s(fishid, k = 19, bs = 're'),
                                       data = data_wels_night_5min,
                                       family = Gamma(link = "log"), select=TRUE,method="REML",
                                       nthreads = 10, cluster = 10, gc.level = 0,
                                       AR.start = startindex, rho = rho_start_value)
```

# <a name="headafinalmodels"></a> Fit random-intercept models with fixed-effects[:page_facing_up:](#headindex)

These models include the fixed effects of the covariates + the random intercept (model_fixed)
``` r
mdl_tench_day_simple_fixed <- bam(formula = det_depth ~
                             s(seasonal_depth, bs = 'ts') +
                             s(amplitude, bs = 'ts') +
                             s(mean_gradient, bs = 'ts') +
                             s(dets_ts5min, bs = 'ts') +
                             s(fishid, k = 19, bs = 're'),
                             data = data_tench_day_5min,
                             family = 'gaussian', select=TRUE,method="REML",
                             nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_tench_day_simple_fixed, plot = TRUE)

#Model with autocorrelation
mdl_tench_day_fixed <- bam(formula = det_depth ~
                                     s(seasonal_depth, bs = 'ts') +
                                     s(amplitude, bs = 'ts') +
                                     s(mean_gradient, bs = 'ts') +
                                     s(dets_ts5min, bs = 'ts') +
                                     s(fishid, k = 19, bs = 're'),
                                     data = data_tench_day_5min,
                                     family = 'gaussian', select=TRUE,method="REML",
                                     nthreads = 10, cluster = 10, gc.level = 0,
                                     AR.start = startindex, rho = rho_start_value)



mdl_tench_night_simple_fixed <- bam(formula = det_depth ~
                                              s(seasonal_depth, bs = 'ts') +
                                              s(amplitude, bs = 'ts') +
                                              s(mean_gradient, bs = 'ts') +
                                              s(dets_ts5min, bs = 'ts') +
                                              s(seasonal_depth, fishid, bs = 're') +
                                              s(amplitude, fishid, bs = 're') +
                                              s(mean_gradient, fishid, bs = 're') +
                                              s(fishid, k = 19, bs = 're'),
                                              data = data_tench_night_5min,
                                              family = Gamma(link = "log"), select=TRUE,method="REML",
                                              nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_tench_night_simple_fixed, plot = TRUE)

#Model with autocorrelation
mdl_tench_night_fixed <- bam(formula = det_depth ~
                             s(seasonal_depth, bs = 'ts') +
                             s(amplitude, bs = 'ts') +
                             s(mean_gradient, bs = 'ts') +
                             s(dets_ts5min, bs = 'ts') +
                             s(seasonal_depth, fishid, bs = 're') +
                             s(amplitude, fishid, bs = 're') +
                             s(mean_gradient, fishid, bs = 're') +
                             s(fishid, k = 19, bs = 're'),
                             data = data_tench_night_5min,
                             family = Gamma(link = "log"), select=TRUE,method="REML",
                             nthreads = 10, cluster = 10, gc.level = 0,
                              AR.start = startindex, rho = rho_start_value)


mdl_pike_day_simple_fixed <- bam(formula = det_depth ~
                             s(seasonal_depth, bs = 'ts') +
                             s(amplitude, bs = 'ts') +
                             s(mean_gradient, bs = 'ts') +
                             s(dets_ts5min, bs = 'ts') +
                             s(fishid, k = 19, bs = 're'),
                             data = data_pike_day_5min,
                             family = 'gaussian', select=TRUE,method="REML",
                             nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_pike_day_simple_fixed, plot = TRUE)

#Model with autocorrelation
mdl_pike_day_fixed <- bam(formula = det_depth ~
                                    s(seasonal_depth, bs = 'ts') +
                                    s(amplitude, bs = 'ts') +
                                    s(mean_gradient, bs = 'ts') +
                                    s(dets_ts5min, bs = 'ts') +
                                    s(fishid, k = 19, bs = 're'),
                                    data = data_pike_day_5min,
                                    family = 'gaussian', select=TRUE,method="REML",
                                    nthreads = 10, cluster = 10, gc.level = 0,
                                    AR.start = startindex, rho = rho_start_value)



mdl_pike_night_simple_fixed <- bam(formula = det_depth ~
                                             s(seasonal_depth, bs = 'ts') +
                                             s(amplitude, bs = 'ts') +
                                             s(mean_gradient, bs = 'ts') +
                                             s(dets_ts5min, bs = 'ts') +
                                             s(seasonal_depth, fishid, bs = 're') +
                                             s(amplitude, fishid, bs = 're') +
                                             s(mean_gradient, fishid, bs = 're') +
                                             s(fishid, k = 12, bs = 're'),
                                             data = data_pike_night_5min,
                                             family = 'gaussian', select=TRUE,method="REML",
                                             nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_pike_night_simple_fixed, plot = TRUE)

#Model with autocorrelation
mdl_pike_night_fixed <- bam(formula = det_depth ~
                                      s(seasonal_depth, bs = 'ts') +
                                      s(amplitude, bs = 'ts') +
                                      s(mean_gradient, bs = 'ts') +
                                      s(dets_ts5min, bs = 'ts') +
                                      s(seasonal_depth, fishid, bs = 're') +
                                      s(amplitude, fishid, bs = 're') +
                                      s(mean_gradient, fishid, bs = 're') +
                                      s(fishid, k = 12, bs = 're'),
                                      data = data_pike_night_5min,
                                      family = 'gaussian', select=TRUE,method="REML",
                                      nthreads = 10, cluster = 10, gc.level = 0,
                                      AR.start = startindex, rho = rho_start_value)


mdl_wels_day_simple_fixed <- bam(formula = det_depth ~
                                           s(seasonal_depth, bs = 'ts') +
                                           s(amplitude, bs = 'ts') +
                                           s(mean_gradient, bs = 'ts') +
                                           s(dets_ts5min, bs = 'ts') +
                                           s(fishid, k = 19, bs = 're'),
                                           data = data_wels_day_5min,
                                           family = Gamma(link = "log"), select=TRUE,method="REML",
                                           nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_wels_day_simple_fixed, plot = TRUE)

#Model with autocorrelation
mdl_wels_day_fixed <- bam(formula = det_depth ~
                                    s(seasonal_depth, bs = 'ts') +
                                    s(amplitude, bs = 'ts') +
                                    s(mean_gradient, bs = 'ts') +
                                    s(dets_ts5min, bs = 'ts') +
                                    s(fishid, k = 19, bs = 're'),
                                    data = data_wels_day_5min,
                                    family = Gamma(link = "log"), select=TRUE,method="REML",
                                    nthreads = 10, cluster = 10, gc.level = 0,
                                    AR.start = startindex, rho = rho_start_value)


mdl_wels_night_simple_fixed <- bam(formula = det_depth ~
                                             s(seasonal_depth, bs = 'ts') +
                                             s(amplitude, bs = 'ts') +
                                             s(mean_gradient, bs = 'ts') +
                                             s(dets_ts5min, bs = 'ts') +
                                             s(fishid, k = 19, bs = 're'),
                                             data = data_wels_night_5min,
                                             family = Gamma(link = "log"), select=TRUE,method="REML",
                                             nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_wels_night_simple_fixed, plot = TRUE)

#Model with autocorrelation
mdl_wels_night_fixed <- bam(formula = det_depth ~
                                      s(seasonal_depth, bs = 'ts') +
                                      s(amplitude, bs = 'ts') +
                                      s(mean_gradient, bs = 'ts') +
                                      s(dets_ts5min, bs = 'ts') +
                                      s(fishid, k = 19, bs = 're'),
                                      data = data_wels_night_5min,
                                      family = Gamma(link = "log"), select=TRUE,method="REML",
                                      nthreads = 10, cluster = 10, gc.level = 0,
                                      AR.start = startindex, rho = rho_start_value)



mdl_rudd_day_simple_fixed <- bam(formula = det_depth ~
                                           s(seasonal_depth, bs = 'ts') +
                                           s(amplitude, bs = 'ts') +
                                           s(mean_gradient, bs = 'ts') +
                                           s(dets_ts5min, bs = 'ts') +
                                           s(fishid, k = 15, bs = 're'),
                                           data = data_rudd_day_5min,
                                           family = 'gaussian', select=TRUE,method="REML",
                                           nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_rudd_day_simple_fixed, plot = TRUE)

#Model with autocorrelation
mdl_rudd_day_fixed <- bam(formula = det_depth ~
                                    s(seasonal_depth, bs = 'ts') +
                                    s(amplitude, bs = 'ts') +
                                    s(mean_gradient, bs = 'ts') +
                                    s(dets_ts5min, bs = 'ts') +
                                    s(fishid, k = 15, bs = 're'),
                                    data = data_rudd_day_5min,
                                    family = 'gaussian', select=TRUE,method="REML",
                                    nthreads = 10, cluster = 10, gc.level = 0,
                                    AR.start = startindex, rho = rho_start_value)


mdl_rudd_night_simple_fixed <- bam(formula = det_depth ~
                                             s(seasonal_depth, bs = 'ts') +
                                             s(amplitude, bs = 'ts') +
                                             s(mean_gradient, bs = 'ts') +
                                             s(dets_ts5min, bs = 'ts') +
                                             s(fishid, k = 15, bs = 're'),
                                             data = data_rudd_night_5min,
                                             family = gaussian(link="log"), select=TRUE,method="REML",
                                             nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_rudd_night_simple_fixed, plot = TRUE)

#Model with autocorrelation
mdl_rudd_night_fixed <- bam(formula = det_depth ~
                                      s(seasonal_depth, bs = 'ts') +
                                      s(amplitude, bs = 'ts') +
                                      s(mean_gradient, bs = 'ts') +
                                      s(dets_ts5min, bs = 'ts') +
                                      s(fishid, k = 15, bs = 're'),
                                      data = data_rudd_night_5min,
                                      family = gaussian(link="log"), select=TRUE,method="REML",
                                      nthreads = 10, cluster = 10, gc.level = 0,
                                      AR.start = startindex, rho = rho_start_value)
```

# <a name="headafinalmodels"></a> Fit random-intercept-random-slopes models [:page_facing_up:](#headindex)

These are final models including fixed effects of the covariates + random intercepts and random slopes for each covariate (model_BRNs).

``` r

mdl_tench_day_simple_BRNs <- bam(formula = det_depth ~
                                           s(seasonal_depth, bs = 'ts') +
                                           s(amplitude, bs = 'ts') +
                                           s(mean_gradient, bs = 'ts') +
                                           s(dets_ts5min, bs = 'ts') +
                                           s(seasonal_depth, fishid, bs = 're') +
                                           s(amplitude, fishid, bs = 're') +
                                           s(mean_gradient, fishid, bs = 're') +
                                           s(fishid, k = 19, bs = 're'),
                                           data = data_tench_day_5min,
                                           family = 'gaussian', select=TRUE,method="REML",
                                           nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_tench_day_simple_BRNs, plot = TRUE)

#Model with autocorrelation
mdl_tench_day_BRNs <- bam(formula = det_depth ~
                                    s(seasonal_depth, bs = 'ts') +
                                    s(amplitude, bs = 'ts') +
                                    s(mean_gradient, bs = 'ts') +
                                    s(dets_ts5min, bs = 'ts') +
                                    s(seasonal_depth, fishid, bs = 're') +
                                    s(amplitude, fishid, bs = 're') +
                                    s(mean_gradient, fishid, bs = 're') +
                                    s(fishid, k = 19, bs = 're'),
                                    data = data_tench_day_5min,
                                    family = 'gaussian', select=TRUE,method="REML",
                                    nthreads = 10, cluster = 10, gc.level = 0,
                                    AR.start = startindex, rho = rho_start_value)


mdl_tench_night_simple_BRNs <- bam(formula = det_depth ~
                                             s(seasonal_depth, bs = 'ts') +
                                             s(amplitude, bs = 'ts') +
                                             s(mean_gradient, bs = 'ts') +
                                             s(dets_ts5min, bs = 'ts') +
                                             s(seasonal_depth, fishid, bs = 're') +
                                             s(amplitude, fishid, bs = 're') +
                                             s(mean_gradient, fishid, bs = 're') +
                                             s(fishid, k = 19, bs = 're'),
                                             data = data_tench_night_5min,
                                             family = Gamma(link = "log"), select=TRUE,method="REML",
                                             nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_tench_night_simple_BRNs, plot = TRUE)

#Model with autocorrelation
mdl_tench_night_BRNs <- bam(formula = det_depth ~
                                      s(seasonal_depth, bs = 'ts') +
                                      s(amplitude, bs = 'ts') +
                                      s(mean_gradient, bs = 'ts') +
                                      s(dets_ts5min, bs = 'ts') +
                                      s(seasonal_depth, fishid, bs = 're') +
                                      s(amplitude, fishid, bs = 're') +
                                      s(mean_gradient, fishid, bs = 're') +
                                      s(fishid, k = 19, bs = 're'),
                                      data = data_tench_night_5min,
                                      family = Gamma(link = "log"), select=TRUE,method="REML",
                                      nthreads = 10, cluster = 10, gc.level = 0,
                                      AR.start = startindex, rho = rho_start_value)


mdl_pike_day_simple_BRNs <- bam(formula = det_depth ~
                                          s(seasonal_depth, bs = 'ts') +
                                          s(amplitude, bs = 'ts') +
                                          s(mean_gradient, bs = 'ts') +
                                          s(dets_ts5min, bs = 'ts') +
                                          s(seasonal_depth, fishid, bs = 're') +
                                          s(amplitude, fishid, bs = 're') +
                                          s(mean_gradient, fishid, bs = 're') +
                                          s(fishid, k = 12, bs = 're'),
                                          data = data_pike_day_5min,
                                          family = 'gaussian', select=TRUE,method="REML",
                                          nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_pike_day_simple_BRNs, plot = TRUE)

#Model with autocorrelation
mdl_pike_day_BRNs <- bam(formula = det_depth ~
                                   s(seasonal_depth, bs = 'ts') +
                                   s(amplitude, bs = 'ts') +
                                   s(mean_gradient, bs = 'ts') +
                                   s(dets_ts5min, bs = 'ts') +
                                   s(seasonal_depth, fishid, bs = 're') +
                                   s(amplitude, fishid, bs = 're') +
                                   s(mean_gradient, fishid, bs = 're') +
                                   s(fishid, k = 12, bs = 're'),
                                   data = data_pike_day_5min,
                                   family = 'gaussian', select=TRUE,method="REML",
                                   nthreads = 10, cluster = 10, gc.level = 0,
                                   AR.start = startindex, rho = rho_start_value)


mdl_pike_night_simple_BRNs <- bam(formula = det_depth ~
                                            s(seasonal_depth, bs = 'ts') +
                                            s(amplitude, bs = 'ts') +
                                            s(mean_gradient, bs = 'ts') +
                                            s(dets_ts5min, bs = 'ts') +
                                            s(seasonal_depth, fishid, bs = 're') +
                                            s(amplitude, fishid, bs = 're') +
                                            s(mean_gradient, fishid, bs = 're') +
                                            s(fishid, k = 12, bs = 're'),
                                            data = data_pike_night_5min,
                                            family = 'gaussian', select=TRUE,method="REML",
                                            nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_pike_night_simple_BRNs, plot = TRUE)

#Model with autocorrelation
mdl_pike_night_BRNs <- bam(formula = det_depth ~
                                     s(seasonal_depth, bs = 'ts') +
                                     s(amplitude, bs = 'ts') +
                                     s(mean_gradient, bs = 'ts') +
                                     s(dets_ts5min, bs = 'ts') +
                                     s(seasonal_depth, fishid, bs = 're') +
                                     s(amplitude, fishid, bs = 're') +
                                     s(mean_gradient, fishid, bs = 're') +
                                     s(fishid, k = 12, bs = 're'),
                                     data = data_pike_night_5min,
                                     family = 'gaussian', select=TRUE,method="REML",
                                     nthreads = 10, cluster = 10, gc.level = 0,
                                     AR.start = startindex, rho = rho_start_value)



mdl_wels_day_simple_BRNs <- bam(formula = det_depth ~
                                          s(seasonal_depth, bs = 'ts') +
                                          s(amplitude, bs = 'ts') +
                                          s(mean_gradient, bs = 'ts') +
                                          s(dets_ts5min, bs = 'ts') +
                                          s(seasonal_depth, fishid, bs = 're') +
                                          s(amplitude, fishid, bs = 're') +
                                          s(mean_gradient, fishid, bs = 're') +
                                          s(fishid, k = 12, bs = 're'),
                                          data = data_wels_day_5min,
                                          family = Gamma(link = "log"), select=TRUE,method="REML",
                                          nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_wels_day_simple_BRNs, plot = TRUE)

#Model with autocorrelation
mdl_wels_day_BRNs <- bam(formula = det_depth ~
                                   s(seasonal_depth, bs = 'ts') +
                                   s(amplitude, bs = 'ts') +
                                   s(mean_gradient, bs = 'ts') +
                                   s(dets_ts5min, bs = 'ts') +
                                   s(seasonal_depth, fishid, bs = 're') +
                                   s(amplitude, fishid, bs = 're') +
                                   s(mean_gradient, fishid, bs = 're') +
                                   s(fishid, k = 12, bs = 're'),
                                   data = data_wels_day_5min,
                                   family = Gamma(link = "log"), select=TRUE,method="REML",
                                   nthreads = 10, cluster = 10, gc.level = 0,
                                   AR.start = startindex, rho = rho_start_value)



mdl_wels_night_simple_BRNs <- bam(formula = det_depth ~
                                            s(seasonal_depth, bs = 'ts') +
                                            s(amplitude, bs = 'ts') +
                                            s(mean_gradient, bs = 'ts') +
                                            s(dets_ts5min, bs = 'ts') +
                                            s(seasonal_depth, fishid, bs = 're') +
                                            s(amplitude, fishid, bs = 're') +
                                            s(mean_gradient, fishid, bs = 're') +
                                            s(fishid, k = 12, bs = 're'),
                                            data = data_wels_night_5min,
                                            family = Gamma(link = "log"), select=TRUE,method="REML",
                                            nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_wels_night_simple_BRNs, plot = TRUE)

#Model with autocorrelation
mdl_wels_night_BRNs <- bam(formula = det_depth ~
                                     s(seasonal_depth, bs = 'ts') +
                                     s(amplitude, bs = 'ts') +
                                     s(mean_gradient, bs = 'ts') +
                                     s(dets_ts5min, bs = 'ts') +
                                     s(seasonal_depth, fishid, bs = 're') +
                                     s(amplitude, fishid, bs = 're') +
                                     s(mean_gradient, fishid, bs = 're') +
                                     s(fishid, k = 12, bs = 're'),
                                     data = data_wels_night_5min,
                                     family = Gamma(link = "log"), select=TRUE,method="REML",
                                     nthreads = 10, cluster = 10, gc.level = 0,
                                     AR.start = startindex, rho = rho_start_value)



mdl_rudd_day_simple_BRNs <- bam(formula = det_depth ~
                                          s(seasonal_depth, bs = 'ts') +
                                          s(amplitude, bs = 'ts') +
                                          s(mean_gradient, bs = 'ts') +
                                          s(dets_ts5min, bs = 'ts') +
                                          s(seasonal_depth, fishid, bs = 're') +
                                          s(amplitude, fishid, bs = 're') +
                                          s(mean_gradient, fishid, bs = 're') +
                                          s(fishid, k = 15, bs = 're'),
                                          data = data_rudd_day_5min,
                                          family = 'gaussian', select=TRUE,method="REML",
                                          nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_rudd_day_simple_BRNs, plot = TRUE)

#Model with autocorrelation
mdl_rudd_day_BRNs <- bam(formula = det_depth ~
                                   s(seasonal_depth, bs = 'ts') +
                                   s(amplitude, bs = 'ts') +
                                   s(mean_gradient, bs = 'ts') +
                                   s(dets_ts5min, bs = 'ts') +
                                   s(seasonal_depth, fishid, bs = 're') +
                                   s(amplitude, fishid, bs = 're') +
                                   s(mean_gradient, fishid, bs = 're') +
                                   s(fishid, k = 15, bs = 're'),
                                   data = data_rudd_day_5min,
                                   family = 'gaussian', select=TRUE,method="REML",
                                   nthreads = 10, cluster = 10, gc.level = 0,
                                   AR.start = startindex, rho = rho_start_value)



mdl_rudd_night_simple_BRNs <- bam(formula = det_depth ~
                                            s(seasonal_depth, bs = 'ts') +
                                            s(amplitude, bs = 'ts') +
                                            s(mean_gradient, bs = 'ts') +
                                            s(dets_ts5min, bs = 'ts') +
                                            s(seasonal_depth, fishid, bs = 're') +
                                            s(amplitude, fishid, bs = 're') +
                                            s(mean_gradient, fishid, bs = 're') +
                                            s(fishid, k = 15, bs = 're'),
                                            data = data_rudd_night_5min,
                                            family=gaussian(link="log"), select=TRUE,method="REML",
                                            nthreads = 10, cluster = 10, gc.level = 0)

#Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_rudd_night_simple_BRNs, plot = TRUE)

#Model with autocorrelation
mdl_rudd_night_BRNs <- bam(formula = det_depth ~
                                     s(seasonal_depth, bs = 'ts') +
                                     s(amplitude, bs = 'ts') +
                                     s(mean_gradient, bs = 'ts') +
                                     s(dets_ts5min, bs = 'ts') +
                                     s(seasonal_depth, fishid, bs = 're') +
                                     s(amplitude, fishid, bs = 're') +
                                     s(mean_gradient, fishid, bs = 're') +
                                     s(fishid, k = 15, bs = 're'),
                                     data = data_rudd_night_5min,
                                     family=gaussian(link="log"), select=TRUE,method="REML",
                                     nthreads = 10, cluster = 10, gc.level = 0,
                                     AR.start = startindex, rho = rho_start_value)
```



# <a name="headafinalmodels"></a>Final models[:page_facing_up:](#headindex)

# <a name="headSpeciesPike"></a>_Pike_[:page_facing_up:](#headindex)

## <a name="headSpeciesPike"></a>![Image name](/outputs/icons/pike_body_1.png) [:sunny:](#head1)


``` r
m1s_pike_day_20k_high_random <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
                                                s(amplitude, bs="cr", k=10) +
                                                s(mean_gradient, bs="cr", k=10) +
                                                s(seasonal_depth, fishid, bs="fs", k=20, m=1) +
                                                s(amplitude, fishid, bs="fs", k=20, m=1) +
                                                s(mean_gradient, fishid, bs="fs", k=20, m=1) +
                                                s(dets_ts, bs="cr", k=10) +
                                                s(fishid, dets_ts, bs="fs", k=20, m=1),
                                                data = data_pike_day,
                                                family = 'gaussian', discrete=TRUE,
                                                nthreads=20, cluster=20, gc.level=0)

rho_start_value <- start_value_rho(m1s_pike_day_20k_high_random, plot = TRUE)

m1_pike_day_20k_high_random <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
                                               s(amplitude, bs="cr", k=10) +
                                               s(mean_gradient, bs="cr", k=10) +
                                               s(seasonal_depth, fishid, bs="fs", k=20, m=1) +
                                               s(amplitude, fishid, bs="fs", k=20, m=1) +
                                               s(mean_gradient, fishid, bs="fs", k=20, m=1) +
                                               s(dets_ts, bs="cr", k=10) +
                                               s(fishid, dets_ts, bs="fs", k=20, m=1),
                                               data = data_pike_day,
                                               family = 'gaussian', discrete=TRUE,
                                               nthreads=20, cluster=20, gc.level=0,
                                               AR.start = startindex, rho = rho_start_value)
```
