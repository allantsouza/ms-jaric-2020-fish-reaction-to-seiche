# <a name="headindex">Influence of internal seiche dynamics on vertical distribution of fish
## <a name="headindex">Seiche MS v5
---

## [Data](#headdata)

[![Image name](/outputs/icons/pike_index_1.png)](#headSpeciesPike)
[![Image name](/outputs/icons/wels_index_1.png)](#headSpeciesWels)
[![Image name](/outputs/icons/rudd_index_1.png)](#headSpeciesRudd)
[![Image name](/outputs/icons/tench_index_1.png)](#headSpeciesTench)

## [Final models table ](#headafinaltable)

## [Random-intercept-only models (model_fishid)](#headamodel_fishid)

## [Random-intercept models with fixed-effects (model_fixed)](#headamodel_fixed)

## [Random-intercept-random-slopes models (model_BRNs)](#headamodel_BRNs)

## [Variance components and Repeatability](#headavars)

## [Tables](#headatables)

## [Plotting reaction norms to seiche and individual variability](#headaplots)

## [References](#headreferences)

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

## <a name="headamodel_fishid"></a> Random-intercept-only models (model_fishid) [:page_facing_up:](#headindex)

These individual-level models serve to calculate the raw phenotypic variance without controlling for fixed-effects of covariates
``` r
# Model without autocorrelation
mdl_tench_day_simple_fishid <- bam(formula = det_depth ~
                                             s(fishid, k = 19, bs = 're'),
                                             data = data_tench_day,
                                             family = 'gaussian', select=TRUE,method="REML",
                                             nthreads = 10, cluster = 10, gc.level = 0)

# Assessing the starting rho value
rho_start_value <- start_value_rho(mdl_tench_day_simple_fishid, plot = TRUE)

# Model with autocorrelation
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

rho_start_value <- start_value_rho(mdl_tench_night_simple_fishid, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_pike_day_simple_fishid, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_pike_night_simple_fishid, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_wels_day_simple_fishid, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_wels_night_simple_fishid, plot = TRUE)

mdl_wels_night_fishid <- bam(formula = det_depth ~
                                       s(fishid, k = 19, bs = 're'),
                                       data = data_wels_night_5min,
                                       family = Gamma(link = "log"), select=TRUE,method="REML",
                                       nthreads = 10, cluster = 10, gc.level = 0,
                                       AR.start = startindex, rho = rho_start_value)
```

## <a name="headamodel_fixed"></a> Random-intercept models with fixed-effects (model_fixed)[:page_facing_up:](#headindex)

These models include the fixed effects of the covariates + the random intercept
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

rho_start_value <- start_value_rho(mdl_tench_day_simple_fixed, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_tench_night_simple_fixed, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_pike_day_simple_fixed, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_pike_night_simple_fixed, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_wels_day_simple_fixed, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_wels_night_simple_fixed, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_rudd_day_simple_fixed, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_rudd_night_simple_fixed, plot = TRUE)

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

## <a name="headamodel_BRNs"></a> Random-intercept-random-slopes models (model_BRNs)[:page_facing_up:](#headindex)

These are final models including fixed effects of the covariates + random intercepts and random slopes for each covariate

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

rho_start_value <- start_value_rho(mdl_tench_day_simple_BRNs, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_tench_night_simple_BRNs, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_pike_day_simple_BRNs, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_pike_night_simple_BRNs, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_wels_day_simple_BRNs, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_wels_night_simple_BRNs, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_rudd_day_simple_BRNs, plot = TRUE)

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

rho_start_value <- start_value_rho(mdl_rudd_night_simple_BRNs, plot = TRUE)

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

## <a name="headavars"></a>Variance components and Repeatability[:page_facing_up:](#headindex)

Create separate lists containing all fitted models
``` r
model_fishid <- list(mdl_tench_day_fishid,
                     mdl_tench_night_fishid,
                     mdl_pike_day_fishid,
                     mdl_pike_night_fishid,
                     mdl_wels_day_fishid,
                     mdl_wels_night_fishid,
                     mdl_rudd_day_fishid,
                     mdl_rudd_night_fishid)

model_fixed <- list(mdl_tench_day_fixed,
                    mdl_tench_night_fixed,
                    mdl_pike_day_fixed,
                    mdl_pike_night_fixed,
                    mdl_wels_day_fixed,
                    mdl_wels_night_fixed,
                    mdl_rudd_day_fixed,
                    mdl_rudd_night_fixed)

model_BRNs <- list(mdl_tench_day_BRNs,
                   mdl_tench_night_BRNs,
                   mdl_pike_day_BRNs,
                   mdl_pike_night_BRNs,
                   mdl_wels_day_BRNs,
                   mdl_wels_night_BRNs,
                   mdl_rudd_day_BRNs,
                   mdl_rudd_night_BRNs)
```

Calculate variance components on the original scale

``` r
for(model in model_BRNs){
        vars <- function(model){
        var_rsc <- c(model$reml.scale/(model$sp[1]/model$smooth[[1]]$S.scale),
                     model$reml.scale/(model$sp[2]/model$smooth[[2]]$S.scale),
                     model$reml.scale/(model$sp[3]/model$smooth[[3]]$S.scale),
                     model$reml.scale/(model$sp[4]/model$smooth[[4]]$S.scale),
                     model$reml.scale/(model$sp[5]/model$smooth[[5]]$S.scale),
                     model$reml.scale/(model$sp[6]/model$smooth[[6]]$S.scale),
                     model$reml.scale/(model$sp[7]/model$smooth[[7]]$S.scale),
                     model$reml.scale/(model$sp[8]/model$smooth[[8]]$S.scale),
                     scale = model$reml.scale)
          vars <- data.frame(var = var_rsc, gam.vcomp(model))
                  print(vars)
}}

```

Calculate repeatabilities
``` r
# Create a list for each species and diel period

names <- c("Tench-Day", "Tench-Night", "Pike-Day", "Pike-Night", "Wels-Day",
    "Wels-Night", "Rudd-Day", "Rudd-Night")


# Function for calculating unadjusted and adjusted R (± 95% CI) using variance components

for(model in model_fishid){
unadj.R <- function(model){
 vars <- variance_comp(model, rescale = TRUE, coverage = 0.95)
 unadj.R <- (vars$variance)[1]/((vars$variance)[1] + (vars$variance)[2])
 unadj.R_lower_ci <- (vars$lower_ci)[1]^2/((vars$lower_ci)[1]^2 + (vars$lower_ci)[2]^2)
 unadj.R_upper_ci <- (vars$upper_ci)[1]^2/((vars$upper_ci)[1]^2 + (vars$upper_ci)[2]^2)
    print(c(unadj.R_lower_ci,unadj.R,unadj.R_upper_ci))

}}

for(model in model_fixed){
adj.fix.R<- function(model){
 vars <- variance_comp(model, rescale = TRUE, coverage = 0.95)
 adj.fix.R <-(vars$variance)[5]/((vars$variance)[5] + (vars$variance)[6])
 adj.fix.R_lower_ci <- (vars$lower_ci)[5]^2/((vars$lower_ci)[6]^2 + (vars$lower_ci)[6]^2)
 adj.fix.R_upper_ci <- (vars$upper_ci)[5]^2/((vars$upper_ci)[6]^2 + (vars$upper_ci)[6]^2)
    print(c(adj.fix.R_lower_ci,adj.fix.R,adj.fix.R_upper_ci))
}}

for(model in model_BRNs){
adj.BRNs.R <- function(model){
 vars <- variance_comp(model, rescale = TRUE, coverage = 0.95)
 adj.BRNs.R <- (vars$variance)[8]/((vars$variance)[8] + (vars$variance)[9])
 adj.BRNs.R_lower_ci <-(vars$lower_ci)[8]^2/((vars$lower_ci)[8]^2 + (vars$lower_ci)[9]^2)
 adj.BRNs.R_upper_ci <-(vars$upper_ci)[8]^2/((vars$upper_ci)[8]^2 + (vars$upper_ci)[9]^2)
    print(c(adj.BRNs.R_lower_ci,adj.BRNs.R, adj.BRNs.R_upper_ci))
}}
```

Create separate tables for unadjusted and adjusted R
``` r
table_unadj.R <- lapply(model_fishid, unadj.R)
df_unadj.R <- as.data.frame(data.table::transpose((table_unadj.R)))
colnames(df_unadj.R) <- c("lower CI", "unadj.R", "upper CI")
unadj_R <- zapsmall(df_unadj.R, digits=6)
unadj_R <- cbind(ID = 1:nrow(unadj_R), unadj_R)    # add index applying cbind function
unadj_R$ID <- as.numeric(unadj_R$ID)

table_adj.fix.R <- lapply(model_fix, adj.fix.R)
df_adj.fix.R <- as.data.frame(data.table::transpose((table_adj.fix.R)))
colnames(df_adj.fix.R) <- c("lower CI", "adj.fix.R", "upper CI")
adj_fix_R <- zapsmall(df_adj.fix.R, digits=6)
adj_fix_R <- cbind(ID = 1:nrow(adj_fix_R), adj_fix_R)
adj_fix_R$ID <- as.numeric(adj_fix_R$ID)

table_adj.BRNs.R <- lapply(model_BRNs, adj.BRNs.R)
df_adj.BRNs.R <- as.data.frame(data.table::transpose((table_adj.BRNs.R)))
colnames(df_adj.BRNs.R) <- c("lower CI", "adj.BRNs.R", "upper CI")
adj_BRNs_R <- zapsmall(df_adj.BRNs.R, digits=6)
adj_BRNs_R <- cbind(ID = 1:nrow(adj_BRNs_R), adj_BRNs_R)
adj_BRNs_R$ID <- as.numeric(adj_BRNs_R$ID)

# Merge the three tables

All_R <- merge(unadj_R, adj_fix_R, by =c("ID"),all=TRUE)
All_R <- merge(All_R, adj_BRNs_R, by =c("ID"),all=TRUE)
rownames(All_R) <- names
```

## <a name="headatables"></a>Tables [:page_facing_up:](#headindex)

|    <br>Predictors                   |    <br>Northern Pike            |                                  |    <br>Wels catfish             |                                  |    <br>Tench                    |                                  |    <br>Rudd                     |                                  |
|-------------------------------------|---------------------------------|----------------------------------|---------------------------------|----------------------------------|---------------------------------|----------------------------------|---------------------------------|----------------------------------|
|    <br>                             |    <br>Day<br>   <br>(n=12)     |    <br>Night<br>   <br>(n=12)    |    <br>Day<br>   <br>(n=15)     |    <br>Night<br>   <br>(n=15)    |    <br>Day<br>   <br>(n=19)     |    <br>Night<br>   <br>(n=19)    |    <br>Day<br>   <br>(n=15)     |    <br>Night<br>   <br>(n=15)    |
|    <br>Parametric   coefficients    |    <br>                         |    <br>                          |    <br>                         |    <br>                          |    <br>                         |    <br>                          |    <br>                         |    <br>                          |
|    <br>(Intercept)                  |    <br>6.80                     |    <br>6.25 ***                  |    <br>1.56 ***                 |    <br>1.06 ***                  |    <br>4.08 ***                 |    <br>0.95 ***                  |    <br>2.99 ***                 |    <br>0.52 ***                  |
|    <br>Smooth terms                 |    <br>                         |    <br>                          |    <br>                         |    <br>                          |    <br>                         |    <br>                          |    <br>                         |    <br>                          |
|    <br>s(seasonal depth)            |    <br>7.88                     |    <br>8.34                      |    <br>8.94 ***                 |    <br>8.78 ***                  |    <br>8.54 ***                 |    <br>8.86 ***                  |    <br>8.52 ***                 |    <br>8.92 ***                  |
|    <br>s(amplitude)                 |    <br>7.33 *                   |    <br>7.53                      |    <br>8.67 **                  |    <br>8.30                      |    <br>5.86 ***                 |    <br>7.38                      |    <br>7.92 **                  |    <br>8.70 ***                  |
|    <br>s(mean gradient)             |    <br>4.07                     |    <br>4.80                      |    <br>8.51 **                  |    <br>7.46                      |    <br>6.93 ***                 |    <br>8.28 *                    |    <br>7.24                     |    <br>8.13 ***                  |
|    <br>s (time)                     |    <br>3.63 ***                 |    <br>8.36 ***                  |    <br>8.98 ***                 |    <br>8.95 ***                  |    <br>8.86 ***                 |    <br>8.91 ***                  |    <br>8.85 ***                 |    <br>8.95 ***                  |
|    <br>s(fishID.seasonal depth)     |    <br>10.79 ***                |    <br>10.73 ***                 |    <br>13.99 ***                |    <br>13.96 ***                 |    <br>14.84 ***                |    <br>15.47 ***                 |    <br>5.13 ***<br>   <br>      |    <br>13.25 ***<br>   <br>      |
|    <br>s(fishID.amplitude)          |    <br>9.40 ***<br>   <br>      |    <br>9.93 ***<br>   <br>       |    <br>13.70 ***<br>   <br>     |    <br>13.41 ***<br>   <br>      |    <br>9.32 **<br>   <br>       |    <br>16.37 ***<br>   <br>      |    <br>9.55 ***<br>   <br>      |    <br>12.15 ***                 |
|    <br>s(fishID.mean gradient)      |    <br>8.66 ***<br>   <br>      |    <br>9.53 ***<br>   <br>       |    <br>13.41 ***<br>   <br>     |    <br>13.04 ***<br>   <br>      |    <br>11.90 ***<br>   <br>     |    <br>16.63 ***<br>   <br>      |    <br>9.96 ***<br>   <br>      |    <br>12.78 ***<br>   <br>      |
|    <br>s(fishID)                    |    <br>10.79 ***<br>   <br>     |    <br>10.81 ***<br>   <br>      |    <br>13.99 ***<br>   <br>     |    <br>13.99 ***<br>   <br>      |    <br>15.99 ***<br>   <br>     |    <br>16.64 ***                 |    <br>12.59 ***<br>   <br>     |    <br>13.32 ***<br>   <br>      |
|    <br>Scale parameter (ϕ)          |    <br>3.86                     |    <br>4.63                      |    <br>0.09                     |    <br>0.32                      |    <br>2.14                     |    <br>0.27                      |    <br>0.96                     |    <br>0.69                      |
|    <br>Adjusted R2                  |    <br>0.54                     |    <br>0.56                      |    <br>0.59                     |    <br>0.56                      |    <br>0.39                     |    <br>0.19                      |    <br>0.35                     |    <br>0.33                      |
|    <br>Deviance explained   (%)     |    <br>53.8                     |    <br>56.5                      |    <br>50.4                     |    <br>48.3                      |    <br>38.7                     |    <br>19.1                      |    <br>35                       |    <br>33                        |

**Table 2** Generalized additive mixed models (GAMM) analyzing trends in depth movement of fish according to physical properties of thermocline (seasonal depth, amplitude and mean gradient) and over time. Numbers refer to estimates (parametric coefficients), estimated effective degrees of freedom (EDF) for main-effects smoothers, and variance components (± 95% confidence intervals, CI) for random-effects smoothers where fishID specifies random factor intercepts and random slopes for each predictor variable. Inter-individual variability is modelled with similar functional responses between fish across context-dependent covariates. In all models, s() are smoothing functions of the covariates computed using penalized regression splines. The maximum number of knots (k-value) used for fitting each model was set to 10 for all splines. EDF indicates the degree of complexity (i.e. wiggliness) of the non-linear relationship between a covariate and the dependent variable, depth. Scale parameter (ϕ), the estimated scale parameter used to compute the variance of i.i.d. Gaussian random effects of the penalized coefficients. Adjusted-R2, measures the goodness-of-fit of the model as the proportion of variance explained. Deviance explained, refers to the proportion of null deviance explained by the model. P-values: * p<0.05, ** p<0.01, *** p<0.001.


|                       |    <br>R0<br>   <br>(unadj)          |    <br>R<br>   <br>(adj fixed)       |    <br>VfishID                          |    <br>Ve                             |    <br>VfishID(sesonal_depth)          |    <br>VfishID(amplitude)               |    <br>VfishID(mean_gradient)           |    <br>R<br>   <br>(adj BRNs)        |
|-----------------------|--------------------------------------|--------------------------------------|-----------------------------------------|---------------------------------------|----------------------------------------|-----------------------------------------|-----------------------------------------|--------------------------------------|
|    <br>Pike-Day       |    <br>0.50<br>   <br>[0.30–0.69]    |    <br>0.50<br>   <br>[0.22–1.15]    |    <br>160.9<br>   <br> [67.4–383.7]    |    <br>3.86<br>   <br> [3.82–3.91]    |    <br>2.20<br>   <br> [0.92–5.24]     |    <br>0.024<br>   <br> [0.018–0.13]    |    <br>0.007<br>   <br> [0.001–0.01]    |    <br>0.98<br>   <br>[0.95–0.99]    |
|    <br>Pike-Night     |    <br>0.62<br>   <br>[0.41–0.79]    |    <br>0.61<br>   <br>[0.34–1.80]    |    <br>4.95<br>   <br>[2.11–11.65]      |    <br>4.65<br>   <br>[4.60–4.72]     |    <br>2.26<br>   <br>[0.95–5.4]       |    <br>0.024<br>   <br>[0.01–0.06]      |    <br>0.007<br>   <br>[0.003–0.02]     |    <br>0.52<br>   <br>[0.32–0.71]    |
|    <br>Wels-Day       |    <br>0.31<br>   <br>[0.18–0.49]    |    <br>0.36<br>   <br>[0.14–0.59]    |    <br>0.06<br>   <br> [0.03–0.12]      |    <br>0.09<br>   <br> [0.09–0.10]    |    <br>0.02<br>   <br> [0.011–0.05]    |    <br>0.001<br>   <br>[0.001–0.002]    |    <br>0.001<br>   <br> [0–0.001]       |    <br>0.40<br>   <br>[0.24–0.58]    |
|    <br>Wels-Night     |    <br>0.32<br>   <br>[0.18–0.49]    |    <br>0.35<br>   <br>[0.13–0.57]    |    <br>0.19<br>   <br>[0.09–0.4]        |    <br>0.32<br>   <br>[0.32–0.33]     |    <br>0.03<br>   <br>[0.014–0.06]     |    <br>0.002<br>   <br>[0.001–0.01]     |    <br>0.001<br>   <br>[0.001–0.003]    |    <br>0.37<br>   <br>[0.22–0.55]    |
|    <br>Tench-Day      |    <br>0.18<br>   <br>[0.10–0.31]    |    <br>0.16<br>   <br>[0.09–0.28]    |    <br>0.49<br>   <br> [0.245–0.99]     |    <br>2.14<br>   <br>[2.12–2.16]     |    <br>0.28<br>   <br>[0.136–0.59]     |    <br>0.003<br>   <br>[0.001–0.01]     |    <br>0.007<br>   <br>[0.003–0.02]     |    <br>0.19<br>   <br>[0.10–0.32]    |
|    <br>Tench-Night    |    <br>0.11<br>   <br>[0.06–0.19]    |    <br>0.15<br>   <br>[0.05–0.17]    |    <br>0.06<br>   <br> [0.03–0.14]      |    <br>0.27<br>   <br>[0.27–0.28]     |    <br>0.03<br>   <br>[0.01–0.06]      |    <br>0.01<br>   <br>[0.004–0.03]      |    <br>0.005<br>   <br>[0.002–0.01]     |    <br>0.19<br>   <br>[0.10–0.33]    |
|    <br>Rudd-Day       |    <br>0.02<br>   <br>[0.01–0.04]    |    <br>0.11<br>   <br>[0.03–0.13]    |    <br>0.11<br>   <br> [0.05–0.26]      |    <br>0.96<br>   <br>[0.95–0.97]     |    <br>0.007<br>   <br>[0.001–0.05]    |    <br>0.002<br>   <br>[0.001–0.005]    |    <br>0.002<br>   <br>[0.001–0.01]     |    <br>0.11<br>   <br>[0.05–0.21]    |
|    <br>Rudd-Night     |    <br>0.03<br>   <br>[0.02–0.07]    |    <br>0.04<br>   <br>[0.01–0.04]    |    <br>0.07<br>   <br> [0.03–0.16]      |    <br>0.69<br>   <br>[0.68–0.70]     |    <br>0.06<br>   <br>[0.03–0.13]      |    <br>0.001<br>   <br>[0.001–0.003]    |    <br>0.002<br>   <br>[0.001–0.004]    |    <br>0.10<br>   <br>[0.05–0.18]    |

**Table 3** Repeatability in individual behavior across contexts from thermocline –dependent variables. Data show unadjusted and adjusted R values with 95% confidence intervals (CI) computed from GAMs model_fishid, model_fixed, and model_BRNs. Variance components are also shown only for model_BRNs. The variances of the random effects for each covariate are calculated by dividing the estimated scale parameter (ϕ) by the smoothness parameter (λ) from the fitted GAMM. Repeatability was calculated using a non-linear approximation to the definition by Nakagawa and Schielzeth (2010) as the relationship between the inter-individual variance (VfishID), given by the smoothing random intercept s(fishID) (i.e., τ00), and the total phenotypic variance (Vtotal = VfishID + Ve), where Ve is the residual variance of the random effects (or intra-individual variability) and equals to ϕ (see SM for more details on variance components and R estimation).



## <a name="headaplots"></a>Plotting reaction norms to seiche and individual variability [:page_facing_up:](#headindex)

:chart: **Figure X**
``` r
tiff("Fig.X.tiff", width = 465, height = 265, units='mm', res = 300)

par(omi=rep(1.0, 4), mar=c(0,3,3,0), mfrow=c(2,4))

plot_smooth(mdl_pike_day_5min_t, view = "amplitude", h0 = c(7.372877), v0 = 0, eegAxis = TRUE, hide.label = TRUE, legend_plot_all = FALSE, plot_all = "fishid", xaxt="n",  ylim=c(0,15), xlim=c(-2,2), xlab = "", ylab = "", main = "", cex.axis = 1.8, cex.lab = 1.8, cex.main =2.5, col = cbPalette[1:20])
axis(1, col.axis = "transparent", lwd = 1)
text(0.5, 12, "zero seiche", cex = 1)
box()
mtexti("Northern pike", 3, cex = 2.8, font = 2)
plot_smooth(mdl_rudd_day_5min_t, view = "amplitude", h0 = c(3.05802), v0 = 0, eegAxis = TRUE, hide.label = TRUE, legend_plot_all = FALSE, plot_all = "fishid",  xaxt="n",  ylim=c(0,4), xlim=c(-2,2), xlab = "", ylab = "", main = " ", cex.axis = 1.8, cex.lab = 1.8, cex.main =2.5, col = cbPalette[1:20])
box()
axis(1, col.axis = "transparent", lwd = 1)
mtexti("Rudd", 3, cex = 2.8, font = 2)
plot_smooth(mdl_tench_day_5min_t, view = "amplitude", h0 = c(4.01068), v0 = 0, eegAxis = TRUE, hide.label = TRUE, legend_plot_all = FALSE, plot_all = "fishid", xaxt="n" ,ylim=c(0,8), xlim=c(-2,2), xlab = "", ylab = "", main = "", cex.axis = 1.8, cex.lab = 1.8, cex.main =2.5, col = cbPalette[1:20])
box()
axis(1, col.axis = "transparent", lwd = 1)
mtexti("Tench", 3, cex = 2.8, font = 2)
plot_smooth(mdl_wels_day_5min_t, view = "amplitude", h0 = c(5.08355), v0 = 0, transform = exp, eegAxis = TRUE, hide.label = TRUE, legend_plot_all = FALSE,  xaxt="n", ylim=c(0,8), xlim=c(-2,2), plot_all = "fishid",  xlab = "", ylab = "", main = "", cex.axis = 2, cex.lab = 2, cex.main =2.5, col = cbPalette[1:20])
box()
axis(1, col.axis = "transparent", lwd = 1)
mtexti("Wels catfish", 3, cex = 2.8, font = 2)
mtexti("Day", 4, cex = 2.8, font = 2)
plot_smooth(mdl_pike_night_5min_t, view = "amplitude", h0 = c(6.613485), v0 = 0, eegAxis = TRUE, hide.label = TRUE, legend_plot_all = FALSE, plot_all = "fishid", ylim=c(0,15), xlim=c(-2,2), xlab = "", ylab = "", main = "", cex.axis = 1.8, cex.lab = 1.8, cex.main =2.5, col = cbPalette[1:20])
box()
plot_smooth(mdl_rudd_night_5min_t, view = "amplitude", h0 = c(1.818448), v0 = 0, transform = exp, eegAxis = TRUE, hide.label = TRUE, legend_plot_all = FALSE, ylim=c(0,4), xlim=c(-2,2), plot_all = "fishid", xlab = "", ylab = "", main = "", cex.axis = 1.8, cex.lab = 1.8, cex.main =2.5, col = cbPalette[1:20])
box()
plot_smooth(mdl_tench_night_5min_t, view = "amplitude", h0 = c(2.692947), v0 = 0, transform = exp, eegAxis = TRUE, hide.label = TRUE, legend_plot_all = FALSE, ylim=c(0,8), xlim=c(-2,2), plot_all = "fishid", xlab = "", ylab = "", main = "", cex.axis = 1.8, cex.lab = 1.8, cex.main =2.5, col = cbPalette[1:20])
box()
plot_smooth(mdl_wels_night_5min_t, view = "amplitude", h0 = c(3.341507), v0 = 0, transform = exp, eegAxis = TRUE, hide.label = TRUE, legend_plot_all = FALSE, ylim=c(0,8), xlim=c(-2,2), plot_all = "fishid",xlab = "", ylab = "", main = "", cex.axis = 1.8, cex.lab = 1.8, cex.main =2.5, col = cbPalette[1:20])
box()
mtexti("Night", 4, cex = 2.8, font = 2)
mtext("amplitude",side=1,line=4,outer=TRUE,cex=1.8)
mtext("Depth (m)",side=2,line=2,outer=TRUE,cex=1.8,las=0)

dev.off()
```
![Seiche_models](/outputs/plots/Fig.X.png "Seiche_models")

**Fig.X** Plots of fish response to internal seiche as a function of gradual thermocline amplitude (i.e., difference of seasonal and actual thermocline depth). Depth was measured repeatedly every one minute (5-min averaged) between June 16 and October 15 in 12 adult northern pike (Esox lucius), 15 rudd (Scardinius erythropthalmus), 19 tench (Tinca tinca) and 15 wels catfish (Silurus glanis) during daytime and at night. Each line corresponds to an individual trajectory across increasing seiche amplitude (left to right) with zero representing periods without seiche (dotted line). Mean depth is denoted by a vertical line. Non-linear trends were modelled as reaction norms, indicative of behavioral plasticity, by including random intercept and random slope effects in the model. The shaded area represents the 95% CI for the fitted curve indicating the range of action of an individual or variability in a reaction norm.


:chart: **Figure Y**
``` r
# Create a common dataset by species-diel_period

all_sp_diel <- droplevels(detections[!detections$species == 'roach',])   # drop unused level
all_sp_diel$diel_period<-as.factor(all_sp_diel$diel_period)

all_sp_diel = all_sp_diel %>%
                             filter(is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))

# Calculate 5-min average and convert into new df

all_sp_diel$dets_ts5min <- floor_date(all_sp_diel$dets_ts, "5 mins")
all_sp_diel_5min <- all_sp_diel %>% group_by(fishid, species, diel_period, dets_ts5min, startindex) %>%
                                              summarise_at(vars(det_depth,seasonal_depth,amplitude,mean_gradient),mean)

all_sp_diel_5min<-as.data.table(all_sp_diel_5min)

# Define amplitude intervals by creating a new categorical varaible amplitude_cat (finally not used)

all_sp_diel_5min[amplitude <= -1.5, amplitude_cat := "-1.5"]
all_sp_diel_5min[amplitude > -1.5 & amplitude <= -1.25, amplitude_cat := "-1.25"]
all_sp_diel_5min[amplitude > -1 & amplitude <= -0.75, amplitude_cat := "-1"]
all_sp_diel_5min[amplitude > -0.75 & amplitude <= -0.50, amplitude_cat := "-0.75"]
all_sp_diel_5min[amplitude > -0.50 & amplitude <= -0.25, amplitude_cat := "-0.50"]
all_sp_diel_5min[amplitude > -0.25 & amplitude <= 0, amplitude_cat := "-0.25"]
all_sp_diel_5min[amplitude > 0 & amplitude <= 0.25, amplitude_cat := "0"]
all_sp_diel_5min[amplitude > 0.25 & amplitude <= 0.5, amplitude_cat := "0.25"]
all_sp_diel_5min[amplitude > 0.5 & amplitude <= 0.75, amplitude_cat := "0.5"]
all_sp_diel_5min[amplitude > 0.75 & amplitude <= 1, amplitude_cat := "0.75"]
all_sp_diel_5min[amplitude > 1 & amplitude <= 1.25, amplitude_cat := "1"]
all_sp_diel_5min[amplitude > 1.25 & amplitude <= 1.5, amplitude_cat := "1.25"]
all_sp_diel_5min[amplitude > 1.5, amplitude_cat := "1.5"]

all_sp_diel_5min <- all_sp_diel_5min %>% group_by(species, diel_period, fishid, amplitude_cat) %>%
                                         summarise_at(vars(det_depth, amplitude), mean)

all_sp_diel_5min<-as.data.frame(all_sp_diel_5min)
all_sp_diel_5min$amplitude_cat<-as.integer(all_sp_diel_5min$amplitude_cat)
all_sp_diel_5min<-all_sp_diel_5min[complete.cases(all_sp_diel_5min), ]
all_sp_diel_5min$fishid = reorder(all_sp_diel_5min$fishid, all_sp_diel_5min$det_depth, FUN = mean)  # increasing order of fish mean depth

# calculate mean depth for each species-diel_period

library(tidyverse)
data_mean <- all_sp_diel_5min %>%
    group_by(species, diel_period) %>%
    summarise(depth = mean(det_depth))

# change names levels (species-diel_period) and define labels

label_species <- as_labeller(c('pike'="Northern pike",
                           'rudd'="Rudd",
                           'tench'="Tench",
                           'wels'="Wels catfish"))

label_diel_period <- as_labeller(c('day'="Day",'night'="Night"))


all_labels<-data.frame(fishid=c(3,3.5,4.3,3.6,3,3.5,4.3,3.6), det_depth=c(10.5,10.5,10.5,10.5,10.5,10.5,10.5,10.5), species = c('pike','rudd','tench','wels','pike','rudd','tench','wels'),
                       diel_period = c('day','day','day','day','night','night','night','night'), label= c("R = 0.98", "R = 0.11", "R = 0.10", "R = 0.40", "R = 0.52", "R = 0.10", "R = 0.19", "R = 0.37"))


mid <- mean(all_sp_diel_5min$amplitude)
P<-ggplot(all_sp_diel_5min, aes(fishid, det_depth)) +
             geom_point(aes(color = amplitude), size = 2) +
             geom_line(linetype = "dotted")+
             geom_line(aes(color=amplitude))+
             geom_hline(data = data_mean, aes(yintercept = depth),color="black", linetype="dashed") +
             facet_grid(diel_period~species,scales="free_x",labeller = labeller(species = label_species, diel_period = label_diel_period))+
             geom_text(data = all_labels,label=all_labels$label, size=6, family="serif", fontface="italic")  +
             theme_bw() +
             theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank()) +
             scale_color_gradient(low = "blue", high = "red",limits = c(-2,2)) +
             scale_y_reverse(limits = c(11, 0)) +
             xlab("Fish ID") + ylab("Depth (m)") + labs(colour = "amplitude") +
             theme(
             plot.title = element_blank(),
             axis.title.x = element_blank(),
             axis.text.x = element_blank(),
             axis.text = element_text( size = 18 ),
             axis.title.y = element_blank(),
             legend.position='right',
             legend.justification='center',
             legend.direction='vertical',
             legend.title = element_text(color = "black", size = 18, face = "plain"),
             legend.text = element_text(color = "black", size = 18),
             legend.key.width=unit(1, "cm"),
             legend.key.height=unit(1, "cm"),
             strip.text = element_text(size = 22, face = "bold"))    # remove face = "bold" for normal grouping labels

y.grob <- textGrob("Depth (m)", gp=gpar(fontface="plain", col="black", fontsize=20), rot=90)
x.grob <- textGrob("", gp=gpar(fontface="plain", col="black", fontsize=20))


tiff("Fig.Y.tiff", width = 465, height = 265, units='mm', res = 300)
p<-grid.arrange(arrangeGrob(P, left = y.grob, bottom = x.grob))
dev.off()
```
![Seiche_models](/outputs/plots/Fig.Y.png "Seiche_models")

**Fig. Y** Inter- and intra-individual variation in depth use across the amplitude of seiche (depicted by color gradient) according to diel periods. Individual fish are arranged in increasing order of mean depth performance (x-axis, left to right). The average depth use of the group is denoted by a dashed horizontal line. A single vertical line denotes variation in a fish’s mean depth across a range of average amplitudes (i.e., individual repeatability) with a longer line indicating higher intra-individual variation. Within a group behavioral differences among fish are shown from left to right, with less overlapping vertical lines indicating overall larger repeatable behavioral differences between the individual fish. Adjusted repeatability (R) was computed from fitted GAMMs after controlling for confounding effects (fixed and random factors) of seasonal depth and mean gradient over time.

:chart: **Figure Z**
``` r

td_Viz <- getViz(mdl_tench_day_BRNs)
tn_Viz <- getViz(mdl_tench_night_BRNs)
pd_Viz <- getViz(mdl_pike_day_BRNs)
pn_Viz <- getViz(mdl_pike_night_BRNs)
wd_Viz <- getViz(mdl_wels_day_BRNs)
wn_Viz <- getViz(mdl_wels_night_BRNs)
rd_Viz <- getViz(mdl_rudd_day_BRNs)
rn_Viz <- getViz(mdl_rudd_night_BRNs)


tiff("Fig.Z.tiff", width = 465, height = 265, units='mm', res = 300)
p<-gridPrint(plot(sm(pd_Viz, 2)) +  theme(
             plot.title = element_text(color = "black", size = 22, face = "bold"),
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text = element_text( size = 14)) +  l_ciLine(colour = 1) + l_ciPoly() + l_fitLine() + scale_y_reverse() + labs(title = "Northern pike", y = "s(amplitude)", x = NULL),
          plot(sm(rd_Viz, 2)) +  theme(
             plot.title = element_text(color = "black", size = 22, face = "bold"),
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank()) +  l_ciLine(colour = 1) + l_ciPoly() + l_fitLine() + scale_y_reverse() + labs(title = "Rudd", x = NULL, y = ""),
          plot(sm(td_Viz, 2)) +  theme(
             plot.title = element_text(color = "black", size = 22, face = "bold"),
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank()) + l_ciLine(colour = 1) + l_ciPoly() + l_fitLine() + scale_y_reverse() + labs(title = "Tench",  x = NULL, y = NULL),
          plot(sm(wd_Viz, 2)) +  theme(
             plot.title = element_text(color = "black", size = 22, face = "bold"),
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank()) +  l_ciLine(colour = 1) + l_ciPoly() + l_fitLine() + scale_y_reverse() + labs(title = "Wels catfish",  x = NULL, y = ""),
          plot(sm(pn_Viz, 2)) +  theme(
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             axis.text = element_text( size = 14)) + l_ciLine(colour = 1) + l_ciPoly() + l_fitLine() + scale_y_reverse() + labs(title = NULL),
          plot(sm(rn_Viz, 2)) +  theme(
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.text = element_text( size = 14)) + l_ciLine(colour = 1) + l_ciPoly() + l_fitLine() + scale_y_reverse() + labs(title = NULL,  y = ""),
          plot(sm(tn_Viz, 2)) +  theme(
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.text = element_text( size = 14)) +  l_ciLine(colour = 1) + l_ciPoly() + l_fitLine() + scale_y_reverse() + labs(title = NULL,  y = ""),
          plot(sm(wn_Viz, 2)) +  theme(
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.text = element_text( size = 14)) +  l_ciLine(colour = 1) + l_ciPoly() + l_fitLine() + scale_y_reverse() + labs(title = NULL,  y = ""),  ncol = 4)

y.grob <- textGrob("s(amplitude)", gp=gpar(fontface="plain", col="black", fontsize=20), rot=90)
x.grob <- textGrob("amplitude", gp=gpar(fontface="plain", col="black", fontsize=20))
right.grob<-textGrob(expression(bold("Day                                   Night")), rot = -90, gp = gpar(col = "black", fontsize = 22))
p<-grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob, right= right.grob))
dev.off()
```
![Seiche_models](/outputs/plots/Fig.Z.png "Seiche_models")

**Fig.Z** Partial effects of the smooth functions of the amplitude of seiche on fish vertical movement of northern pike, rudd, tench and wels catfish according to diel periods. The slopes represent average effects for the amplitude of thermocline based on GAMs fit to 5-min average fish position data over a whole time period. Splines for other predictors are omitted.

