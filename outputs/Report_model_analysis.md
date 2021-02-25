# <a name="headindex">Influence of internal seiche dynamics on vertical distribution of fish
## <a name="headindex">Report of model analysis
---

# [Data](#headdata)

# [Global model formula](#headglobalmodel)

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
string with 14-24 unique values (depending on data subset)
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

```

# <a name="headglobalmodel"></a>Global model formula [:page_facing_up:](#headindex)

The global model formula incorporates an **AR(1)** autoregressive residual term. This is found to be suitable for the different data subsets when compared to others autocorrelation structures. First, we fit the model without AR(1) component in order to estimate the **rho-value** and then re-fit the model with the AR1 structure.<br />
Included is a squared first derivative penalty term (`m=1`) for the _random factor-smooth interactions_ to correct uncertainty around the mean of the _main-effects_ smoothers (with second derivative penalty, `m=2` by default) to reduce concurvity between both smoothers.<br />
Adding the argument `discrete=TRUE` substantially reduces computation time by creating bins of discrete data for each variable before fitting the model.<br />

The global model formula can be expressed as follows:

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{E\left&space;(&space;y_{i}&space;\right&space;)}=&space;\beta_{0}&space;&plus;&space;f_{1}\left&space;(&space;X_{1i}&space;\right&space;)&space;&plus;&space;f_{2}\left&space;(&space;X_{1i}&space;\right&space;)_{id_{i}}&space;&plus;&space;f_{3}\left&space;(&space;X_{2i}&space;\right&space;)&space;&plus;&space;f_{4}\left&space;(&space;X_{2i}&space;\right&space;)_{id_{i}}&space;&plus;&space;f_{5}\left&space;(&space;X_{3i}&space;\right&space;)&space;&plus;&space;f_{6}\left&space;(&space;X_{3i}&space;\right&space;)_{id_{i}}&space;&plus;&space;f_{7}\left&space;(&space;t_{i}&space;\right&space;)&space;&plus;&space;f_{8}\left&space;(&space;t_{i}&space;\right&space;)_{id_{i}}&space;&plus;&space;\varepsilon&space;_{i}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{E\left&space;(&space;y_{i}&space;\right&space;)}=&space;\beta_{0}&space;&plus;&space;f_{1}\left&space;(&space;X_{1i}&space;\right&space;)&space;&plus;&space;f_{2}\left&space;(&space;X_{1i}&space;\right&space;)_{id_{i}}&space;&plus;&space;f_{3}\left&space;(&space;X_{2i}&space;\right&space;)&space;&plus;&space;f_{4}\left&space;(&space;X_{2i}&space;\right&space;)_{id_{i}}&space;&plus;&space;f_{5}\left&space;(&space;X_{3i}&space;\right&space;)&space;&plus;&space;f_{6}\left&space;(&space;X_{3i}&space;\right&space;)_{id_{i}}&space;&plus;&space;f_{7}\left&space;(&space;t_{i}&space;\right&space;)&space;&plus;&space;f_{8}\left&space;(&space;t_{i}&space;\right&space;)_{id_{i}}&space;&plus;&space;\varepsilon&space;_{i}" title="\boldsymbol{E\left ( y_{i} \right )}= \beta_{0} + f_{1}\left ( X_{1i} \right ) + f_{2}\left ( X_{1i} \right )_{id_{i}} + f_{3}\left ( X_{2i} \right ) + f_{4}\left ( X_{2i} \right )_{id_{i}} + f_{5}\left ( X_{3i} \right ) + f_{6}\left ( X_{3i} \right )_{id_{i}} + f_{7}\left ( t_{i} \right ) + f_{8}\left ( t_{i} \right )_{id_{i}} + \varepsilon _{i}" /></a>

where:

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{E\left&space;(&space;y_{i}&space;\right&space;)}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{E\left&space;(&space;y_{i}&space;\right&space;)}" title="\boldsymbol{E\left ( y_{i} \right )}" /></a>
is the expected value of depth for a fish <a href="https://www.codecogs.com/eqnedit.php?latex=_{id_{i}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?_{id_{i}}" title="_{id_{i}}" /></a>
at time <a href="https://www.codecogs.com/eqnedit.php?latex=t_{i}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?t_{i}" title="t_{i}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\beta_{0}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\beta_{0}" title="\beta_{0}" /></a> is the is the average value of the response (model intercept).

<a href="https://www.codecogs.com/eqnedit.php?latex=X_{1i},&space;X_{2i},&space;X_{3i}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?X_{1i},&space;X_{2i},&space;X_{3i}" title="X_{1i}, X_{2i}, X_{3i}" /></a> are main-effects of the continuous covariates (i.e., seasonal depth, mean gradient and amplitude).

<a href="https://www.codecogs.com/eqnedit.php?latex=t_{i}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?t_{i}" title="t_{i}" /></a> is the average time effects trend.

<a href="https://www.codecogs.com/eqnedit.php?latex=f_{1-8}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?f_{1-8}" title="f_{1-8}" /></a> are smoothing functions of the covariates <a href="https://www.codecogs.com/eqnedit.php?latex=X_{1i},&space;X_{2i},&space;X_{3i},&space;t" target="_blank"><img src="https://latex.codecogs.com/svg.latex?X_{1i},&space;X_{2i},&space;X_{3i},&space;t" title="X_{1i}, X_{2i}, X_{3i}, t" /></a> for a fish <a href="https://www.codecogs.com/eqnedit.php?latex=_{id_{i}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?_{id_{i}}" title="_{id_{i}}" /></a>.

<a href="https://www.codecogs.com/eqnedit.php?latex=_{id_{i}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?_{id_{i}}" title="_{id_{i}}" /></a> is the _Fish ID_ group-level factor.

<a href="https://www.codecogs.com/eqnedit.php?latex=\varepsilon&space;_{i}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\varepsilon&space;_{i}" title="\varepsilon _{i}" /></a> is the random error term including the random smooths and an AR(1) autoregressive residual term.


### <a name="head6"></a>:balance_scale: _Model 1: Main-effects smoothers + Group-level smoothers with same wiggliness (Random factor-smoother interaction)_

Depth is modelled as a function of eight smoother terms where we find main-effects functions for the covariates plus fish-specific deviations around that main functions. Thus, change in fish depth is non-linear across the range of each covariate (seasonal depth, amplitude, mean gradient) and over time, all of which are modelled as _factor-smooth interactions_ using argument `bs"fs"`.<br />
With this specification, gam treats random effects as smooths allowing a separate smooth for each level of _fishid_, with the same smoothing parameter for all smooths (see `?mgcv::smooth.construct.fs.smooth.spec`):
- 4 *main-effects* smoothers (with cubic splines) for all observations of each covariate (individual effects have a common penalty term).
- 0 group-specific *random-effect smoothers* (without fish-specific intercepts).
- 4 group-level *random factor-smoother interactions* (i.e., for all _Fish ID_ with `bs="fs"`):
  - A penalty to shrink these smoothers toward zero.
  - A common smoothing parameter for all _Fish ID_ smoothers with same wiggliness (˜<a href="https://www.codecogs.com/eqnedit.php?latex=\simeq" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\simeq" title="\simeq" /></a> functional responses for all _Fish ID_ )
  - Different shapes of the smooth terms (~inter-individual variation in responses, t2).
``` r
m1 <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
		                  s(amplitude, bs="cr", k=10) +
                      s(mean_gradient, bs="cr", k=10) +
                      s(dets_ts, bs="cr", k=10) +
                      s(seasonal_depth, fishid, bs="fs", k=10, m=1) +
                      s(amplitude, fishid, bs="fs", k=10, m=1) +
                      s(mean_gradient, fishid, bs="fs", k=10, m=1) +
                      s(dets_ts, fishid, bs="fs", k=10, m=1) +
                      data=data_wels_day, method="REML", family="gaussian",
                      discrete=TRUE, nthreads=10, cluster=10, gc.level=0,
                      AR.start=startindex, rho=rho_start_value)
```
This is equivalent to:
``` r
m1 <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
		                  s(amplitude, bs="cr", k=10) +
                      s(mean_gradient, bs="cr", k=10) +
                      s(dets_ts, bs="cr", k=10) +
                      t2(seasonal_depth, fishid, bs=c("tp","re"), k=10, m=1, full=TRUE) +
                      t2(amplitude, fishid, bs=c("tp","re"), k=10, m=1, full=TRUE) +
                      t2(mean_gradient, fishid, bs=c("tp","re"), k=10, m=1, full=TRUE) +
                      t2(dets_ts, fishid, bs=c("tp","re"), k=10, m=1, full=TRUE),
                      data=data_wels_day, method="REML", family="gaussian",
                      discrete=TRUE, nthreads=10, cluster=10, gc.level=0,
                      AR.start=startindex, rho=rho_start_value)
```
(but see **Variable selection and parameters estimation of the main and random smooth effects**)

:information_source: _Note that the default spline used is `"tp"` unless otherwise specified. In the models formuia I specify cyclic splines (`"cr"`) for the _main-effects smoothers_ while the _random smoothers_ use the default `"tp"` that is equivalent to specifying `xt="tp"` in the first formula._


# <a name="headafinalmodels"></a>Final models[:page_facing_up:](#headindex)

Based on the global structure (**Model 1**), model selection went through the following steps:

1. Choose the error distribution that fits each data subset best, based on maximum goodness-of-fit estimation, AIC and inspection of residual patterns using QQ plots.
2. Splines selection: try different splines functions (e.g., `"cr"` vs. `"tp"`).<br />
3.  Try different number of basis dimensions (**k-value**). Shift k-value from 10 to 20 to 50 to n:<br />
  - The model fitted uses group-level _random factor-smooth interactions_ where the k-value is determined using the `gam.check()` function.<br />
  - After fitting an initial model with k=10 (default), `gam.check()` is used to see whether more wiggliness is necessary, i.e., whether the smooths use up all of df (edf <a href="https://www.codecogs.com/eqnedit.php?latex=\approx" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\approx" title="\approx" /></a> k´) or not (edf < k`). <br />
  - If `gam.check()` suggests that more wiggliness is necessary (edf <a href="https://www.codecogs.com/eqnedit.php?latex=\approx" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\approx" title="\approx" /></a> k`), this procedure is repeated again using models with increased k-values allowing for the splines to adapt to the wiggly and non-wiggly parts.<br />
  - The k-value can be different for both _main-effects_ smoothers and _random factor-smooth interactions_ (mixed k: more or less wiggliness for each term).<br />
4. Compare models using residuals QQ plots, AIC and <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a>. However, the AIC-based method should be treated with care when applied to models including an AR1 structure or *random factor-smoother interactions* as it is the case. Instead, performing Chi-square test with the functions `anova()` (_mgcv_ package) and `compareML()` (_itsadug_ package) is generally preferred and more accurate (mostly for non-normal error models).<br />
5. Where appropiate (linear relationships), compute a simplified equivalent model by dropping the smoothing functions. Compare the full-smoothing and the simplified models (I omit this step in this report).<br />
6. On the selected model, extract and round the edf of the smooth terms (estimated with penalization).<br />
7. Re-fit the model with fixed edf (un-penalized) for the _main-effects smoothers_ using `fx=TRUE`, excluding the _random group-level smoothers_.<br />
8. Compare models with and without fixed edf.


# <a name="headSpeciesPike"></a>_Pike_[:page_facing_up:](#headindex)

## <a name="headSpeciesPike"></a>![Image name](/outputs/icons/pike_body_1.png) [:sunny:](#head1)

For this specific data set, given that AIC model selection is a questionable method here as outlined in **Final models**, we computed models using either Gaussian or an alternative scaled T-distribution (`"scat"`) with cubic splines (`"cr"`) for the _random factor-smoother interactions_. Overall, the data fitted well to a Gaussian distribution while the latter distribution didn´t improve the model fit.<br />
Models k10, k20 and k50 with equal or mixed k-values (higher or lower for the random terms) were compared through Chi-square test using the function _compareML()_. In all cases, the models presented edf/k´ ratios far from ideal, however, the residuals patterns distribution for the smooth terms of those models were relatively acceptable.<br />
The mixed k50 model (higher for the random terms) was preferred but given that more data are necessary to estimate the larger number of coefficients from the increased complexity, it couldn´t be computed with fixed df. We finally kept the next closest in preference based on Chi-square test and distribution of residual patterns, with a <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a> indicating very good performance.<br />

:balance_scale: **Model**

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
:1234: **Table**

<table border=0>
<caption align="bottom"> SUMMARY OF m1_pike_day_20k_high_random </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 3.18 </td> <td align="right"> 1.11 </td> <td align="right"> 2.88 </td> <th align="right"> 0.004 </th> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 1 </td> <td align="right"> 1 </td> <td align="right"> 5.03 </td> <th align="right"> 0.025 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 1.03 </td> <td align="right"> 1.03 </td> <td align="right"> 0.04 </td> <td align="right"> 0.874 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 1.01 </td> <td align="right"> 1.01 </td> <td align="right"> 0.96 </td> <td align="right"> 0.326 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 1 </td> <td align="right"> 1 </td> <td align="right"> 1.77 </td> <td align="right"> 0.184 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 199.62 </td> <td align="right"> 236 </td> <td align="right"> 19.71 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 126.04 </td> <td align="right"> 230 </td> <td align="right"> 5.42 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 104.12 </td> <td align="right"> 239 </td> <td align="right"> 1.55 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 214.3 </td> <td align="right"> 239 </td> <td align="right"> 30.86 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.80 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 236337.7 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 80.2% </td> </tr>
   <a name=tab.gam></a>
</table>

:chart: **Plots**

![Seiche_models](/outputs/plots/m1_pike_day_20k_high_random.png "Seiche_models")

Re-fit the model with fixed degreees of freedom

:balance_scale: **Final model**
``` r
rho_start_value <- start_value_rho(m1s_pike_day_20k_high_random, plot=TRUE)

f.df <- round(summary(m1_pike_day_20k_high_random)$edf)+1 # get edf per smooth
f.df <- pmax(f.df,3)                                      # minimum basis dimension is 3

m1_pike_day_20k_high_random_final <- bam(det_depth ~ s(seasonal_depth, k=f.df[1], fx=TRUE) +
                                                     s(amplitude, k=f.df[2], fx=TRUE)+
                                                     s(mean_gradient, k=f.df[3], fx=TRUE) +
                                                     s(seasonal_depth, fishid, bs="fs", k=f.df[4], m=1) +
                                                     s(amplitude, fishid, bs="fs", k=f.df[5], m=1) +
                                                     s(mean_gradient, fishid, bs="fs", k=f.df[6], m=1) +
                                                     s(dets_ts, k=f.df[7], fx=TRUE)+
                                                     s(fishid, dets_ts, bs="fs", k=f.df[8], m=1),
                                                     data = data_pike_day,
                                                     family = 'gaussian', discrete=TRUE,
                                                     nthreads=10, cluster=10, gc.level=0,
                                                     AR.start = startindex, rho = rho_start_value)

```
:1234: **Table**

<table border=2>
<caption align="bottom"> SUMMARY OF m1_pike_day_20k_high_random_final </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 5.64 </td> <td align="right"> 2.55 </td> <td align="right"> 2.21 </td> <th align="right"> 0.027 </th> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 2 </td> <td align="right"> 2 </td> <td align="right"> 0.05 </td> <td align="right"> 0.952 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 2 </td> <td align="right"> 2 </td> <td align="right"> 0.64 </td> <td align="right"> 0.525 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 2 </td> <td align="right"> 2 </td> <td align="right"> 1.85 </th> <td align="right"> 0.158 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 2 </td> <td align="right"> 2 </td> <td align="right"> 0.24 </td> <td align="right"> 0.786 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 904.9 </td> <td align="right"> 2064 </td> <td align="right"> 1.73 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 378.1 </td> <td align="right"> 1154 </td> <td align="right"> 0.81 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 139.1 </td> <td align="right"> 1133 </td> <td align="right"> 0.2 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 1787 </td> <td align="right"> 2309 </td> <td align="right"> 12.89 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.92 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 204822.7 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 91.7% </td> </tr>
   <a name=tab.gam></a>
</table>

:chart: **Plots**

![Seiche_models](/outputs/plots/m1_pike_day_20k_high_random_final.png "Seiche_models")


## <a name="headSpeciesPike"></a>![Image name](/outputs/icons/pike_body_1.png) [:waning_crescent_moon:](#head1)

The Gaussian distribution, besides the poor residuals patterns, was an optimal fit based on maximum goodness-of-fit estimation. Models were re-fitted with a scaled T-distribution (`"scat"`) and compared using Chi-square test that showed an improvement in residuals patterns.<br />
Of all models fitted, the one with k=10 was the most robust in terms of deviance and residuals patterns and was always preferred to k20 and k50 models either with equal or mixed k-values (higher or lower for the random terms).<br />
After fixing the df, the k10 model was approximated by an un-penalized mixed k20 model (higher for the random terms) with cubic splines (`"cr"`) instead of  (`"tp"`) basis functions, and higher <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a> (+8) but the Chi-square test determined that the simpler k10 model was preferred.<br />

:balance_scale: **Model**

``` r
m1_pike_night_10k_scat <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
                                          s(amplitude, bs="cr", k=10) +
                                          s(mean_gradient, bs="cr", k=10) +
                                          s(seasonal_depth, fishid, bs="fs", k=10, m=1) +
                                          s(amplitude, fishid, bs="fs", k=10, m=1) +
                                          s(mean_gradient, fishid, bs="fs", k=10, m=1) +
                                          s(dets_ts, bs="cr", k=10) +
                                          s(fishid, dets_ts, bs="fs", k=10, m=1),
                                          data=data_pike_night,
                                          family='scat', discrete=TRUE,
                                          nthreads=10, cluster=10, gc.level=0)

rho_start_value <- start_value_rho(m1s_pike_night_10k_scat, plot=TRUE)

m1_pike_night_10k_scat <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
                                          s(amplitude, bs="cr", k=10) +
                                          s(mean_gradient, bs="cr", k=10) +
                                          s(seasonal_depth, fishid, bs="fs", k=10, m=1) +
                                          s(amplitude, fishid, bs="fs", k=10, m=1) +
                                          s(mean_gradient, fishid, bs="fs", k=10, m=1) +
                                          s(dets_ts, bs="cr", k=10) +
                                          s(fishid, dets_ts, bs="fs", k=10, m=1),
                                          data=data_pike_night,
                                          family='scat', discrete=TRUE,
                                          nthreads=10, cluster=10, gc.level=0,
                                          AR.start = startindex, rho = rho_start_value)

```
:1234: **Table**

<table border=0>
<caption align="bottom"> SUMMARY OF m1_pike_night_10k_scat </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 3.81 </td> <td align="right"> 1.96 </td> <td align="right"> 1.94 </td> <td align="right"> 0.052 </td> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 7.58 </td> <td align="right"> 7.88 </td> <td align="right"> 0.24 </td> <td align="right"> 0.985 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 8.41 </td> <td align="right"> 8.61 </td> <td align="right"> 1.46 </td> <td align="right"> 0.103 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 3.14 </td> <td align="right"> 3.52 </td> <td align="right"> 2.22 </td> <td align="right"> 0.143 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 8.37 </td> <td align="right"> 8.59 </td> <td align="right"> 3.37 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 96.67 </td> <td align="right"> 120 </td> <td align="right"> 31.77 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 77.9 </td> <td align="right"> 120 </td> <td align="right"> 7.55 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 79.85 </td> <td align="right"> 120 </td> <td align="right"> 1.02 </th> <td align="right"> 0.926 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 101.86 </td> <td align="right"> 120 </td> <td align="right"> 27.61 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.67 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 140798.3 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 60.8% </td> </tr>
   <a name=tab.gam></a>
</table>

:chart: **Plots**

![Seiche_models](/outputs/plots/m1_pike_night_10k_scat.png "Seiche_models")

Re-fit the model with fixed degreees of freedom

:balance_scale: **Final model**
``` r
rho_start_value <- start_value_rho(m1s_pike_night_10k_scat, plot=TRUE)

f.df <- round(summary(m1_pike_night_10k_scat)$edf)+1 # get edf per smooth
f.df <- pmax(f.df,3)                                 # minimum basis dimension is 3

m1_pike_night_10k_scat_final <- bam(det_depth ~ s(seasonal_depth, k=f.df[1], fx=TRUE) +
                                                s(amplitude, k=f.df[2], fx=TRUE)+
                                                s(mean_gradient, k=f.df[3], fx=TRUE) +
                                                s(seasonal_depth, fishid, bs="fs", k=f.df[4], m=1) +
                                                s(amplitude, fishid, bs="fs", k=f.df[5], m=1) +
                                                s(mean_gradient, fishid, bs="fs", k=f.df[6], m=1) +
                                                s(dets_ts, k=f.df[7], fx=TRUE)+
                                                s(fishid, dets_ts, bs="fs", k=f.df[8], m=1),
                                                data=data_pike_night,
                                                family='scat', discrete=TRUE,
                                                nthreads=10, cluster=10, gc.level=0,
                                                AR.start = startindex, rho = rho_start_value)
```
:1234: **Table**

<table border=2>
<caption align="bottom"> SUMMARY OF m1_pike_night_10k_scat_final </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 4.36 </td> <td align="right"> 2.66 </td> <td align="right"> 1.64 </td> <td align="right"> 0.101 </td> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 8 </td> <td align="right"> 8 </td> <td align="right"> 0.84 </td> <td align="right"> 0.566 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 8 </td> <td align="right"> 8 </td> <td align="right"> 5.36 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 3 </td> <td align="right"> 3 </td> <td align="right"> 0.15 </th> <td align="right"> 0.93 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 8 </td> <td align="right"> 8 </td> <td align="right"> 0.2 </td> <td align="right"> 0.991 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 709 </td> <td align="right"> 1075 </td> <td align="right"> 12.1 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 454.7 </td> <td align="right"> 780 </td> <td align="right"> 3.72 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 276.7 </td> <td align="right"> 942 </td> <td align="right"> 0.7 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 1012.6 </td> <td align="right"> 1181 </td> <td align="right"> 77.51 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.85 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 53694.64 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 75.4% </td> </tr>
   <a name=tab.gam></a>
</table>


:chart: **Plots**

![Seiche_models](/outputs/plots/m1_pike_night_10k_scat_final.png "Seiche_models")


# <a name="headSpeciesWels"></a>_Wels catfish_[:page_facing_up:](#headindex)

## <a name="headSpeciesWels"></a>![Image name](/outputs/icons/wels_body_1.png) [:sunny:](#head1)

Overall, the data fitted better to a scaled T-distribution where random smooths with cubic splines (`"cr"`) seem to do much better than with (`"tp"`) basis functions, based on AIC and the look of the residuals patterns for all smoothers. However, for the same k-value, models fitted with the Gamma distribution showed an improvement in the edf/k´ ratio as well as a 2-units increase in the <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a> value.<br />
The Gamma distribution was preferred to any other distribution based on maximum goodness-of-fit estimation but models with the two distributions were further compared. The Chi-square tests show that models with the Gamma distribution are indeed a better fit to the data. Of these, k20 models either with equal or mixed k-values (higher for the random terms) exhibit a balance between edf/k´ ratio and <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a> while re-fitting with un-penalized smoother terms unlike models with k <a href="https://www.codecogs.com/eqnedit.php?latex=\geq" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\geq" title="\geq" /></a> 40, having more coefficients to estimate than actual data. The latter models had lower <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a> than fixed k20 models and worse edf/k´ ratios as the k-value increased.<br />
Of the k20 models, the one with higher k-value for the random terms was selected in Chi-square test. In this model, three of the _main-effects_ terms showed a linear relationship with the response.<br />
``` r

AIC_table <- AIC(m1_wels_day_gamma_20k_high_global,m1_wels_day_gamma_20k_high_random,
                 m1_wels_day_gamma_20k,m1_wels_day_gamma_40k,m1_wels_day_gamma_10k)%>%
             rownames_to_column(var= "Model")%>%
             mutate(data_source = rep(c("data_wels_day"), each=1))%>%
             group_by(data_source)%>%
             mutate(deltaAIC = AIC - min(AIC))%>%
             ungroup()%>%
             dplyr::select(-data_source)%>%
             mutate_at(.vars = vars(df,AIC, deltaAIC),
                       .funs = funs(round,.args = list(digits=0)))

kable(AIC_table)
```
|Model                             |   df|    AIC| deltaAIC|
|:---------------------------------|----:|------:|--------:|
|m1_wels_day_gamma_20k_high_global |  375| 902357|      930|
|m1_wels_day_gamma_20k_high_random |  631| 901427|        0|
|m1_wels_day_gamma_20k             |  636| 903617|     2190|
|m1_wels_day_gamma_10k             |  368| 903605|     2178|
|m1_wels_day_gamma_40k             | 1205| 903748|     2322|


:balance_scale: **Model**

``` r
m1s_wels_day_gamma_20k_high_random <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
                                                      s(amplitude, bs="cr", k=10) +
                                                      s(mean_gradient,  bs="cr", k=10) +
                                                      s(seasonal_depth, fishid, bs="fs", k=20, m=1) +
                                                      s(amplitude, fishid, bs="fs", k=20, m=1) +
                                                      s(mean_gradient, fishid, bs="fs", k=20, m=1) +
                                                      s(dets_ts, bs="cr", k=10) +
                                                      s(fishid, dets_ts, bs="fs", k=20, m=1),
                                                      data=data_wels_day,
                                                      family=Gamma(link = "log"), discrete=TRUE,
                                                      nthreads=10, cluster=10, gc.level=0)

rho_start_value <- start_value_rho(m1s_wels_day_gamma_20k_high_random, plot=TRUE)

m1_wels_day_gamma_20k_high_random <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
                                                     s(amplitude, bs="cr", k=10) +
                                                     s(mean_gradient,  bs="cr", k=10) +
                                                     s(seasonal_depth, fishid, bs="fs", k=20, m=1) +
                                                     s(amplitude, fishid, bs="fs", k=20, m=1) +
                                                     s(mean_gradient, fishid, bs="fs", k=20, m=1) +
                                                     s(dets_ts, bs="cr", k=10) +
                                                     s(fishid, dets_ts, bs="fs", k=20, m=1),
                                                     data=data_wels_day,
                                                     family=Gamma(link = "log"), discrete=TRUE,
                                                     nthreads=10, cluster=10, gc.level=0,
                                                     AR.start=startindex, rho=rho_start_value)
```
:1234: **Table**

<table border=0>
<caption align="bottom"> SUMMARY OF m1_wels_day_gamma_20k_high_random </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 1.52 </td> <td align="right"> 0.11 </td> <td align="right"> 14.2 </td> <th align="right"> &lt; 0.001 </th> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 1.03 </td> <td align="right"> 1.03 </td> <td align="right"> 1.13 </td> <td align="right"> 0.295 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 1.02 </td> <td align="right"> 1.02 </td> <td align="right"> 1.13 </td> <td align="right"> 0.290 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 3.77 </td> <td align="right"> 4.56 </td> <td align="right"> 4.09 </td> <th align="right"> 0.002 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 1.37 </td> <td align="right"> 1.4 </td> <td align="right"> 0.38 </td> <td align="right"> 0.491 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 214.94 </td> <td align="right"> 300 </td> <td align="right"> 6.2 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 93.42 </td> <td align="right"> 298 </td> <td align="right"> 1.15 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 53.63 </td> <td align="right"> 300 </td> <td align="right"> 0.33 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 256.85 </td> <td align="right"> 300 </td> <td align="right"> 19.39 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.67 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 901427 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 62.3% </td> </tr>
   <a name=tab.gam></a>
</table>

:chart: **Plots**

![Seiche_models](/outputs/plots/m1_wels_day_gamma_20k_high_random.png "Seiche_models")

Re-fit the model with fixed degreees of freedom

:balance_scale: **Final model**
``` r
rho_start_value <- start_value_rho(m1s_wels_day_gamma_20k_high_random, plot=TRUE)

f.df <- round(summary(m1_wels_day_gamma_20k_high_random)$edf)+1 # get edf per smooth
f.df <- pmax(f.df,3)                                            # minimum basis dimension is 3

m1_wels_day_gamma_20k_high_random_final <- bam(det_depth ~ s(seasonal_depth, k=f.df[1], fx=TRUE) +
                                                           s(amplitude, k=f.df[2], fx=TRUE)+
                                                           s(mean_gradient, k=f.df[3], fx=TRUE) +
                                                           s(seasonal_depth, fishid, bs="fs", k=f.df[4], m=1) +
                                                           s(amplitude, fishid, bs="fs", k=f.df[5], m=1) +
                                                           s(mean_gradient, fishid, bs="fs", k=f.df[6], m=1) +
                                                           s(dets_ts, k=f.df[7], fx=TRUE)+
                                                           s(fishid, dets_ts, bs="fs", k=f.df[8], m=1),
                                                           data=data_wels_day,
                                                           family="scat", discrete=TRUE,
                                                           nthreads=10, cluster=10, gc.level=0,
                                                           AR.start=startindex, rho=rho_start_value)
```
:1234: **Table**

<table border=2>
<caption align="bottom"> SUMMARY OF m1_wels_day_gamma_20k_high_random_final </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 1.59 </td> <td align="right"> 0.43 </td> <td align="right"> 3.66 </td> <th align="right"> 0.003 </th> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 2 </td> <td align="right"> 2 </td> <td align="right"> 0.45 </td> <td align="right"> 0.640 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 2 </td> <td align="right"> 2 </td> <td align="right"> 0.69 </td> <td align="right"> 0.5 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 4 </td> <td align="right"> 4 </td> <td align="right"> 5.2 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 2 </td> <td align="right"> 2 </td> <td align="right"> 0.03 </td> <td align="right"> 0.974 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 572.3 </td> <td align="right"> 2895 </td> <td align="right"> 0.51 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 270.68 </td> <td align="right"> 1190 </td> <td align="right"> 0.46 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 63.25 </td> <td align="right"> 808 </td> <td align="right"> 0.11 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 2871.52 </td> <td align="right"> 3744 </td> <td align="right"> 6.77 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.78 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 659663.6 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 78.4% </td> </tr>
   <a name=tab.gam></a>
</table>


:chart: **Plots**

![Seiche_models](/outputs/plots/m1_wels_day_gamma_20k_high_random_final.png "Seiche_models")


## <a name="headSpeciesWels"></a>![Image name](/outputs/icons/wels_body_1.png) [:waning_crescent_moon:](#head1)

This data fitted better with a Gamma distribution and reached an optimal k-value at 20 when higher k-value for the main-effects terms was specified.<br />
Using Chi-square test, the mixed k20 model (higher for the main-effects) was always preferred to default k20 and mixed k20 (higher for the random terms) models and to k50 and k100 models, irrespective of the k-value category (mixed, non-mixed). Thus, k20 models with mixed k-values (higher for the main-effects) was selected against k50 and k100 models with mixed k-values (higher for the main-effects).<br />
Preference for higher wiggliness of the _main-effects smoothers_ indicates that mean functional response of the population tends to be more variable than at the individual-level, hence, fish-specific deviations around that main-effects functions are lower than expected from average population variability.

:balance_scale: **Model**

``` r
m1s_wels_night_gamma_20k_high_main <- bam(det_depth ~ s(seasonal_depth, k=20, bs="cr") +
                                                      s(amplitude, k=20, bs="cr") +
                                                      s(mean_gradient, k=20, bs="cr") +
                                                      s(seasonal_depth, fishid, bs="fs", k=10, m=1) +
                                                      s(amplitude, fishid, bs="fs", k=10, m=1) +
                                                      s(mean_gradient, fishid, bs="fs", k=10, m=1) +
                                                      s(dets_ts, k=20, bs="cr") +
                                                      s(fishid, dets_ts, bs="fs", k=10, m=1),
                                                      data=data_wels_night,
                                                      family=Gamma(link = "log"), discrete=TRUE,
                                                      nthreads=10, cluster=10, gc.level=0)

rho_start_value <- start_value_rho(m1s_wels_night_gamma_20k_high_main, plot=TRUE)

m1_wels_night_gamma_20k_high_main <- bam(det_depth ~ s(seasonal_depth, k=20, bs="cr") +
                                                     s(amplitude, k=20, bs="cr") +
                                                     s(mean_gradient, k=20, bs="cr") +
                                                     s(seasonal_depth, fishid, bs="fs", k=10, m=1) +
                                                     s(amplitude, fishid, bs="fs", k=10, m=1) +
                                                     s(mean_gradient, fishid, bs="fs", k=10, m=1) +
                                                     s(dets_ts, k=20, bs="cr") +
                                                     s(fishid, dets_ts, bs="fs", k=10, m=1),
                                                     data=data_wels_night,
                                                     family=Gamma(link = "log"), discrete=TRUE,
                                                     nthreads=10, cluster=10, gc.level=0,
                                                     AR.start=startindex, rho=rho_start_value)
```

:1234: **Table**

<table border=0>
<caption align="bottom"> SUMMARY OF m1_wels_night_gamma_20k_high_main </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 1.45 </td> <td align="right"> 0.32 </td> <td align="right"> 4.59 </td> <th align="right"> &lt; 0.001 </th> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 16.62 </td> <td align="right"> 17.67 </td> <td align="right"> 4.81 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 13.42 </td> <td align="right"> 15.48 </td> <td align="right"> 3.93 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 1.01 </td> <td align="right"> 1.01 </td> <td align="right"> 0.29 </td> <td align="right"> 0.61 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 17.3 </td> <td align="right"> 18.35 </td> <td align="right"> 7.45 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 120.61 </td> <td align="right"> 150 </td> <td align="right"> 14.61 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 87.83 </td> <td align="right"> 150 </td> <td align="right"> 5.55 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 91.6 </td> <td align="right"> 150 </td> <td align="right"> 3.22 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 116.79 </td> <td align="right"> 150 </td> <td align="right"> 9.87 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.64 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 428213 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 54.6% </td> </tr>
   <a name=tab.gam></a>
</table>

:chart: **Plots**

![Seiche_models](/outputs/plots/m1_wels_night_gamma_20k_high_main.png "Seiche_models")

Re-fit the model with fixed degreees of freedom

:balance_scale: **Final model**
``` r
rho_start_value <- start_value_rho(m1s_wels_night_gamma_20k_high_main, plot = TRUE)

f.df <- round(summary(m1_wels_night_gamma_20k_high_main)$edf)+1 # get edf per smooth
f.df <- pmax(f.df,3) # minimum basis dimension is 3

m1_wels_night_gamma_20k_high_main_final <- bam(det_depth ~ s(seasonal_depth, k=f.df[1], fx=TRUE) +
                                                           s(amplitude, k=f.df[2], fx=TRUE)+
                                                           s(mean_gradient, k=f.df[3], fx=TRUE) +
                                                           s(seasonal_depth, fishid, bs="fs", k=f.df[4], m=1) +
                                                           s(amplitude, fishid, bs="fs", k=f.df[5], m=1) +
                                                           s(mean_gradient, fishid, bs="fs", k=f.df[6], m=1) +
                                                           s(dets_ts, k=f.df[7], fx=TRUE)+
                                                           s(fishid, dets_ts, bs="fs", k=f.df[8], m=1),
                                                           data=data_wels_night,
                                                           family=Gamma(link = "log"), discrete=TRUE,
                                                           nthreads=10, cluster=10, gc.level=0,
                                                           AR.start=startindex, rho =rho_start_value)
```
:1234: **Table**

<table border=2>
<caption align="bottom"> SUMMARY OF m1_wels_night_gamma_20k_high_main_final </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 1.64 </td> <td align="right"> 0.54 </td> <td align="right"> 3.03 </td> <th align="right"> 0.002 </th> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 17 </td> <td align="right"> 17 </td> <td align="right"> 1.81 </td> <th align="right"> 0.022 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 13 </td> <td align="right"> 13 </td> <td align="right"> 1.56 </td> <td align="right"> 0.09 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 2 </td> <td align="right"> 2 </td> <td align="right"> 0.44 </td> <td align="right"> 0.644 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 17 </td> <td align="right"> 17 </td> <td align="right"> 0.4 </td> <td align="right"> 0.986 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 736.2 </td> <td align="right"> 1754 </td> <td align="right"> 2.49 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 297.5 </td> <td align="right"> 1140 </td> <td align="right"> 0.67 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 209.4 </td> <td align="right"> 1310 </td> <td align="right"> 0.35 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 1369.6 </td> <td align="right"> 1770 </td> <td align="right"> 6.07 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.76 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </td> <td align="right"> 330326.1 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 68.1% </td> </tr>
   <a name=tab.gam></a>
</table>


:chart: **Plots**

![Seiche_models](/outputs/plots/m1_wels_night_gamma_20k_high_main_final.png "Seiche_models")


# <a name="headSpeciesTench"></a>_Tench_[:page_facing_up:](#headindex)

## <a name="headSpeciesTench"></a>![Image name](/outputs/icons/tench_body_1.png) [:sunny:](#head1)

Increasing the k-value from k10 models to k100 models increased <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a> but the residual QQ plots and the edf/k´ ratio indicated bad residual patterns in most smooth terms, as denoted by a higher AIC value.<br />
Rather, computing the model with an alternative scaled T-distribution (`"scat"`) with cubic splines (`"cr"`) for the _random factor-smoother interactions_ overall improved the model fit based on the residual patterns and AIC values. In addition, the scat model was also selected in the  Chi-square tests.<br />
On the other hand, the preference for cubic splines in the random smooth terms could indicate that the inter-individual variation follows cyclic patterns across the range of values of the covariates.<br />
Taking the base model with k=10 and increasing to k=20 appeared to slightly improve the model fit based on the lower AIC, higher <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a> and k-index values. However, the Chi-square test always selected the k10 model when compared to k20 models with equal or mixed k-values (higher or lower for the random terms), hence the former was retained.<br />
The re-fitted k10 scat model with fixed df was preferred to either penalized k50 and k100 models, showing better residuals patterns and edf/k´ ratios for all smooth terms.<br />

``` r
AIC_table <- AIC(m1_tench_day_10k,m1_tench_day_scat,m1_tench_day_fs_cr,m1_tench_day,m1_tench_day_fs_cr_scat)%>%
             rownames_to_column(var= "Model")%>%
             mutate(data_source = rep(c("data_tench_day"), each=1))%>%
             group_by(data_source)%>%
             mutate(deltaAIC = AIC - min(AIC))%>%
             ungroup()%>%
             dplyr::select(-data_source)%>%
             mutate_at(.vars = vars(df,AIC, deltaAIC),
                       .funs = funs(round,.args = list(digits=0)))

kable(AIC_table)
```
|Model                   |  df|     AIC| deltaAIC|
|:-----------------------|---:|-------:|--------:|
|m1_tench_day_10k        | 406| 1059123|     9881|
|m1_tench_day_scat       | 416| 1051297|     2055|
|m1_tench_day_fs_cr      | 406| 1059308|    10066|
|m1_tench_day            | 410| 1059108|     9865|
|m1_tench_day_fs_cr_scat | 408| 1049242|        0|


:balance_scale: **Model**

``` r
m1s_tench_day_fs_cr_scat <- bam(det_depth ~ s(seasonal_depth, k=10) +
                                            s(amplitude, k=10) +
                                            s(mean_gradient, k=10) +
                                            s(seasonal_depth, fishid, bs="fs", xt="cr", k=10, m=1) +
                                            s(amplitude, fishid, bs="fs", xt="cr", k=10, m=1) +
                                            s(mean_gradient, fishid, bs="fs", xt="cr", k=10, m=1) +
                                            s(dets_ts, k=10) +
                                            s(fishid, dets_ts, bs="fs", xt="cr", k=10, m=1),
                                            data=data_tench_day,
                                            family="scat", discrete=TRUE,
                                            nthreads=10, cluster=10, gc.level=0)

rho_start_value <- start_value_rho(m1s_tench_day_fs_cr_scat, plot=TRUE)

m1_tench_day_fs_cr_scat <- bam(det_depth ~ s(seasonal_depth, k=10) +
                                           s(amplitude, k=10) +
                                           s(mean_gradient, k=10) +
                                           s(seasonal_depth, fishid, bs="fs", xt="cr", k=10, m=1) +
                                           s(amplitude, fishid, bs="fs", xt="cr", k=10, m=1) +
                                           s(mean_gradient, fishid, bs="fs", xt="cr", k=10, m=1) +
                                           s(dets_ts, k=10) +
                                           s(fishid, dets_ts, bs="fs", xt="cr", k=10, m=1),
                                           data=data_tench_day,
                                           family="scat", discrete=TRUE,
                                           nthreads=10, cluster=10, gc.level=0,
                                           AR.start=startindex, rho=rho_start_value)
```

:1234: **Table**

<table border=0>
<caption align="bottom"> SUMMARY OF m1_tench_day_fs_cr_scat </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 1.04 </td> <td align="right"> 2.08 </td> <td align="right"> 0.5 </td> <td align="right"> 0.617 </td> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 5.58 </td> <td align="right"> 5.94 </td> <td align="right"> 1.19 </td> <td align="right"> 0.245 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 5.3 </td> <td align="right"> 5.93 </td> <td align="right"> 2.07 </td> <td align="right"> 0.055 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 6.45 </td> <td align="right"> 7.36 </td> <td align="right"> 2.61 </td> <th align="right"> 0.009 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 8.8 </td> <td align="right"> 8.84 </td> <td align="right"> 166.96 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 102.1 </td> <td align="right"> 144 </td> <td align="right"> 15.08 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 77.68 </td> <td align="right"> 174 </td> <td align="right"> 2.29 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 64.25 </td> <td align="right"> 174 </td> <td align="right"> 2.33 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 127.72 </td> <td align="right"> 156 </td> <td align="right"> 45.92 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.47 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 1049242 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 43.5% </td> </tr>
   <a name=tab.gam></a>
</table>

:chart: **Plots**

![Seiche_models](/outputs/plots/m1_tench_day_fs_cr_scat.png "Seiche_models")

Re-fit the model with fixed degreees of freedom

:balance_scale: **Final model**
``` r
rho_start_value <- start_value_rho(m1s_tench_day_fs_cr_scat, plot=TRUE)

f.df <- round(summary(m1_tench_day_fs_cr_scat)$edf)+1 # get edf per smooth
f.df <- pmax(f.df,3)                                  # minimum basis dimension is 3

m1_tench_day_fs_cr_scat_final <- bam(det_depth ~ s(seasonal_depth, k=f.df[1], fx=TRUE) +
                                                 s(amplitude, k=f.df[2], fx=TRUE)+
                                                 s(mean_gradient, k=f.df[3], fx=TRUE) +
                                                 s(seasonal_depth, fishid, bs="fs", k=f.df[4], xt="cr", m=1) +
                                                 s(amplitude, fishid, bs="fs", k=f.df[5], xt="cr", m=1) +
                                                 s(mean_gradient, fishid, bs="fs", k=f.df[6], xt="cr", m=1) +
                                                 s(dets_ts, k=f.df[7], fx=TRUE)+
                                                 s(fishid, dets_ts, bs="fs", k=f.df[8], xt="cr", m=1),
                                                 data=data_tench_day,
                                                 family="scat", discrete=TRUE,
                                                 nthreads=10, cluster=10, gc.level=0,
                                                 AR.start=startindex, rho=rho_start_value)
```
:1234: **Table**

<table border=2>
<caption align="bottom"> SUMMARY OF m1_tench_day_fs_cr_scat_final </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 3.55 </td> <td align="right"> 2.16 </td> <td align="right"> 1.64 </td> <td align="right"> 0.101 </td> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 6 </td> <td align="right"> 6 </td> <td align="right"> 0.25 </td> <td align="right"> 0.958 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 5 </td> <td align="right"> 5 </td> <td align="right"> 1.22 </td> <td align="right"> 0.299 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 6 </td> <td align="right"> 6 </td> <td align="right"> 2.5 </td> <th align="right"> 0.02 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 9 </td> <td align="right"> 9 </td> <td align="right"> 0.18 </td> <td align="right"> 0.997 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 766.54 </td> <td align="right"> 1121 </td> <td align="right"> 8.35 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 110.64 </td> <td align="right"> 1110 </td> <td align="right"> 0.5 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 84.55 </td> <td align="right"> 1027 </td> <td align="right"> 0.44 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 1464.94 </td> <td align="right"> 1650 </td> <td align="right"> 24.8 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.63 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 864749.9 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 55.6% </td> </tr>
   <a name=tab.gam></a>
</table>


:chart: **Plots**

![Seiche_models](/outputs/plots/m1_tench_day_fs_cr_scat_final.png "Seiche_models")


## <a name="headSpeciesTench"></a>![Image name](/outputs/icons/tench_body_1.png) [:waning_crescent_moon:](#head1)

Model selection based on knots shifting, AIC and Chi-square test shows that the default k10 model is the best model. Re-fitting k20, k50 and k100 models with unpenalized terms is only possible if _factor-smooth interaction_ terms are not involved, evidencing the default model with lower k is more consistent regarding the estimated edf and p-values.<br />

``` r
AIC_table <- AIC(m1_tench_night_gamma_20k_high_global,m1_tench_night_gamma_20k_high_random,m1_tench_night_gamma_20k,m1_tench_night_gamma_50k_high_random,m1_tench_night_gamma_50k_high_global,
                 m1_tench_night_gamma_100k_high_random,m1_tench_night_gamma_100k_high_global,m1_tench_night_gamma_100k,m1_tench_night_gamma_50k,m1_tench_night_gamma_10k,m1_tench_night_gamma_10k_fs_cr)%>%
                 rownames_to_column(var= "Model")%>%
                 mutate(data_source = rep(c("data_tench_night"), each=1))%>%
                 group_by(data_source)%>%
                 mutate(deltaAIC = AIC - min(AIC))%>%
                 ungroup()%>%
                 dplyr::select(-data_source)%>%
                 mutate_at(.vars = vars(df,AIC, deltaAIC),
                           .funs = funs(round,.args = list(digits=0)))

kable(AIC_table)
```

|Model                                 |   df|    AIC| deltaAIC|
|:-------------------------------------|----:|------:|--------:|
|m1_tench_night_gamma_10k              |  441| 302121|        0|
|m1_tench_night_gamma_10k_fs_cr        | 407|  308486|     6366|
|m1_tench_night_gamma_20k_high_main    |  448| 304677|     2556|
|m1_tench_night_gamma_20k_high_random  |  754| 311882|     9762|
|m1_tench_night_gamma_20k              |  767| 315813|    13692|
|m1_tench_night_gamma_50k              | 1598| 319621|    17500|
|m1_tench_night_gamma_50k_high_random  | 1579| 314057|    11937|
|m1_tench_night_gamma_50k_high_main    |  824| 315181|    13060|
|m1_tench_night_gamma_100k             | 2617| 314945|    12825|
|m1_tench_night_gamma_100k_high_random | 2556| 306562|     4442|
|m1_tench_night_gamma_100k_high_main   |  920| 307955|     5835|


:balance_scale: **Model**

``` r
m1s_tench_night_gamma_10k <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
                                             s(amplitude, bs="cr", k=10) +
                                             s(mean_gradient, bs="cr", k=10) +
                                             s(seasonal_depth, fishid, bs="fs", k=10, m=1) +
                                             s(amplitude, fishid, bs="fs", k=10, m=1) +
                                             s(mean_gradient, fishid, bs="fs", k=10, m=1) +
                                             s(dets_ts, bs="cr", k=10) +
                                             s(fishid, dets_ts, bs="fs", k=10, m=1),
                                             data=data_tench_night,
                                             family=Gamma(link = "log"), discrete=TRUE,
                                             nthreads=10, cluster=10, gc.level=0)

rho_start_value <- start_value_rho(m1s_tench_night_gamma_10k, plot=TRUE)

m1_tench_night_gamma_10k <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
                                            s(amplitude, bs="cr", k=10) +
                                            s(mean_gradient,  bs="cr", k=10) +
                                            s(seasonal_depth, fishid, bs="fs", k=10, m=1) +
                                            s(amplitude, fishid, bs="fs", k=10, m=1) +
                                            s(mean_gradient, fishid, bs="fs", k=10, m=1) +
                                            s(dets_ts, bs="cr", k=10) +
                                            s(fishid, dets_ts, bs="fs", k=10, m=1),
                                            data=data_tench_night,
                                            family=Gamma(link = "log"), discrete=TRUE,
                                            nthreads=10, cluster=10, gc.level=0,
                                            AR.start=startindex, rho=rho_start_value)
```

:1234: **Table**

<table border=0>
<caption align="bottom"> SUMMARY OF m1_tench_night_gamma_10k </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 1.93 </td> <td align="right"> 0.44 </td> <td align="right"> 4.42 </td> <th align="right"> &lt; 0.001 </th> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 7.70 </td> <td align="right"> 7.99 </td> <td align="right"> 3.58 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 4.62 </td> <td align="right"> 5.13 </td> <td align="right"> 1.36 </td> <td align="right"> 0.219 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 6.12 </td> <td align="right"> 6.796 </td> <td align="right"> 4.24 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 8.79 </td> <td align="right"> 8.87 </td> <td align="right"> 28.90 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 99.62 </td> <td align="right"> 176.00 </td> <td align="right"> 8.88 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 100.48 </td> <td align="right"> 187.00 </td> <td align="right"> 4.29 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 89.99 </td> <td align="right"> 190.00 </td> <td align="right"> 2.31 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 117.09 </td> <td align="right"> 169.00 </td> <td align="right"> 11.88 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.34 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 302121 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 34% </td> </tr>
   <a name=tab.gam></a>
</table>

:chart: **Plots**

![Seiche_models](/outputs/plots/m1_tench_night_gamma_10k.png "Seiche_models")

Re-fit the model with fixed degreees of freedom

:balance_scale: **Final model**
``` r
rho_start_value <- start_value_rho(m1s_tench_night_gamma_10k, plot=TRUE)

f.df <- round(summary(m1_tench_night_gamma_10k)$edf)+1  # get edf per smooth
f.df <- pmax(f.df,3)                                    # minimum basis dimension is 3

m1_tench_night_gamma_10k_final <- bam(det_depth ~ s(seasonal_depth, k=f.df[1], fx=TRUE) +
                                                  s(amplitude, k=f.df[2], fx=TRUE) +
                                                  s(mean_gradient, k=f.df[3], fx=TRUE) +
                                                  s(seasonal_depth, fishid, bs="fs", k=f.df[4], m=1) +
                                                  s(amplitude, fishid, bs="fs", k=f.df[5], m=1) +
                                                  s(mean_gradient, fishid, bs="fs", k=f.df[6], m=1) +
                                                  s(dets_ts, k=f.df[7], fx=TRUE) +
                                                  s(fishid, dets_ts, bs="fs", k=f.df[8], m=1),
                                                  data=data_tench_night,
                                                  family=Gamma(link = "log"), discrete=TRUE,
                                                  nthreads=10, cluster=10, gc.level=0,
                                                  AR.start=startindex, rho=rho_start_value)
```
:1234: **Table**

<table border=2>
<caption align="bottom"> SUMMARY OF m1_tench_night_gamma_10k_final </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 0.92 </td> <td align="right"> 0.64 </td> <td align="right"> 1.44 </td> <td align="right"> 0.15 </td> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 8 </td> <td align="right"> 8 </td> <td align="right"> 0.31 </td> <td align="right"> 0.96 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 5 </td> <td align="right"> 5 </td> <td align="right"> 3.06 </td> <th align="right"> 0.01 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 6 </td> <td align="right"> 6 </td> <td align="right"> 6.15 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 9 </td> <td align="right"> 9 </td> <td align="right"> 0.06 </td> <td align="right"> 0.999 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 629.4 </td> <td align="right"> 1245 </td> <td align="right"> 2.61 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 448 </td> <td align="right"> 1471 </td> <td align="right"> 0.78 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 129.5 </td> <td align="right"> 1497 </td> <td align="right"> 0.34 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 1139.9 </td> <td align="right"> 1504 </td> <td align="right"> 6.7 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.60 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 188051 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 59.3% </td> </tr>
   <a name=tab.gam></a>
</table>


:chart: **Plots**

![Seiche_models](/outputs/plots/m1_tench_night_gamma_10k_final.png "Seiche_models")


# <a name="headSpeciesRudd"></a>_Rudd_[:page_facing_up:](#headindex)

## <a name="headSpeciesRudd"></a>![Image name](/outputs/icons/rudd_body_1.png) [:sunny:](#head1)

While k20 models, with either default or mixed k-values (higher for the random terms) are equally plausible based on AIC (deltaAIC = 1), the non-mixed k20 model is preferred using Chi-square test. However, after re-fitting models with fixed df, the mixed k20 model is selected by AIC and Chi-square test, even though both have very similar values of <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a>.<br />
The selected mixed k20 model (higher for the random terms) can be re-fitted using linear parametric terms for _seasonal depth_ and _mean gradient_ (edf <a href="https://www.codecogs.com/eqnedit.php?latex=\approx" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\approx" title="\approx" /></a> 1; note the linearity of both terms in the plots) but the selection based on AIC shows that keeping the smooth terms is preferred (deltaAIC = 58).<br />
On the other hand, a mixed k100 model (higher for the random terms) has a <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a> value slightly lower than for fixed k20 models and cannot be unpenalized as  denoted by the lower proportion of deviance explained.<br />

``` r
AIC_table <- AIC(m1_rudd_day_20k_high_global,m1_rudd_day_20k_high_random,m1_rudd_day_20k,m1_rudd_day_10k)%>%
             rownames_to_column(var= "Model")%>%
             mutate(data_source = rep(c("data_rudd_day"), each=1))%>%
             group_by(data_source)%>%
             mutate(deltaAIC = AIC - min(AIC))%>%
             ungroup()%>%
             dplyr::select(-data_source)%>%
             mutate_at(.vars = vars(df,AIC, deltaAIC),
                       .funs = funs(round,.args = list(digits=0)))

kable(AIC_table)
```
|Model                       |  df|    AIC| deltaAIC|
|:---------------------------|---:|------:|--------:|
|m1_rudd_day_20k_high_global | 327| 218026|     1375|
|m1_rudd_day_20k_high_random | 546| 216652|        1|
|m1_rudd_day_20k             | 551| 216651|        0|
|m1_rudd_day_10k             | 331| 217820|     1169|


``` r
AIC_table <- AIC(m1_rudd_day_20k_high_random,m1_rudd_day_20k,m1_rudd_day_20k_final,m1_rudd_day_20k_high_random_final)%>%
             rownames_to_column(var= "Model")%>%
             mutate(data_source = rep(c("data_rudd_day"), each=1))%>%
             group_by(data_source)%>%
             mutate(deltaAIC = AIC - min(AIC))%>%
             ungroup()%>%
             dplyr::select(-data_source)%>%
             mutate_at(.vars = vars(df,AIC, deltaAIC),
                       .funs = funs(round,.args = list(digits=0)))

kable(AIC_table)
```
|Model                             |   df|    AIC| deltaAIC|
|:---------------------------------|----:|------:|--------:|
|m1_rudd_day_20k_high_random       |  546| 216652|    14937|
|m1_rudd_day_20k                   |  551| 216651|    14935|
|m1_rudd_day_20k_final             | 2237| 202020|      305|
|m1_rudd_day_20k_high_random_final | 2254| 201716|        0|


:balance_scale: **Model**

``` r
m1s_rudd_day_20k_high_random <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
                                                s(amplitude, bs="cr", k=10) +
                                                s(mean_gradient,  bs="cr", k=10) +
                                                s(seasonal_depth, fishid, bs="fs", k=20, m=1) +
                                                s(amplitude, fishid, bs="fs", k=20, m=1) +
                                                s(mean_gradient, fishid, bs="fs", k=20, m=1) +
                                                s(dets_ts, bs="cr", k=10) +
                                                s(fishid, dets_ts, bs="fs", k=20, m=1),
                                                data=data_rudd_day,
                                                family="gaussian", discrete=TRUE,
                                                nthreads=10, cluster=10, gc.level=0)

rho_start_value <- start_value_rho(m1s_rudd_day_20k_high_random, plot=TRUE)

m1_rudd_day_20k_high_random <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
                                               s(amplitude, bs="cr", k=10) +
                                               s(mean_gradient,  bs="cr", k=10) +
                                               s(seasonal_depth, fishid, bs="fs", k=20, m=1) +
                                               s(amplitude, fishid, bs="fs", k=20, m=1) +
                                               s(mean_gradient, fishid, bs="fs", k=20, m=1) +
                                               s(dets_ts, bs="cr", k=10) +
                                               s(fishid, dets_ts, bs="fs", k=20, m=1),
                                               data=data_rudd_day,
                                               family="gaussian", discrete=TRUE,
                                               nthreads=10, cluster=10, gc.level=0,
                                               AR.start=startindex, rho=rho_start_value)
```
:1234: **Table**

<table border=0>
<caption align="bottom"> SUMMARY OF m1_rudd_day_20k_high_random </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 2.09 </td> <td align="right"> 0.69 </td> <td align="right"> 3.03 </td> <th align="right"> 0.002 </th> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 1.41 </td> <td align="right"> 1.46 </td> <td align="right"> 0.14 </td> <td align="right"> 0.879 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 4.89 </td> <td align="right"> 5.35 </td> <td align="right"> 1.29 </td> <td align="right"> 0.346 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 1.01 </td> <td align="right"> 1.02 </td> <td align="right"> 0.23 </td> <td align="right"> 0.64 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 6.77 </td> <td align="right"> 6.9 </td> <td align="right"> 0.79 </td> <td align="right"> 0.566 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 144.61 </td> <td align="right"> 247 </td> <td align="right"> 7.38 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 125.16 </td> <td align="right"> 290 </td> <td align="right"> 1.99 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 83.12 </td> <td align="right"> 286 </td> <td align="right"> 0.72 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 173.43 </td> <td align="right"> 257 </td> <td align="right"> 19.53 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.57 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 216652 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 57.5% </td> </tr>
   <a name=tab.gam></a>
</table>

:chart: **Plots**

![Seiche_models](/outputs/plots/m1_rudd_day_20k_high_random.png "Seiche_models")


Re-fit the model with fixed degreees of freedom

:balance_scale: **Final model**
``` r
rho_start_value <- start_value_rho(m1s_rudd_day_20k_high_random, plot=TRUE)

f.df <- round(summary(m1_rudd_day_20k_high_random)$edf)+1 # get edf per smooth
f.df <- pmax(f.df,3)                                      # minimum basis dimension is 3

m1_rudd_day_20k_high_random_final <- bam(det_depth ~ s(seasonal_depth, k=f.df[1], fx=TRUE) +
                                                     s(amplitude, k=f.df[2], fx=TRUE) +
                                                     s(mean_gradient, k=f.df[3], fx=TRUE) +
                                                     s(seasonal_depth, fishid, bs="fs", k=f.df[4], m=1) +
                                                     s(amplitude, fishid, bs="fs", k=f.df[5], m=1) +
                                                     s(mean_gradient, fishid, bs="fs", k=f.df[6], m=1) +
                                                     s(dets_ts, k=f.df[7], fx=TRUE) +
                                                     s(fishid, dets_ts, bs="fs", k=f.df[8], m=1),
                                                     data=data_rudd_day,
                                                     family="gaussian", discrete=TRUE,
                                                     nthreads=10, cluster=10, gc.level=0,
                                                     AR.start=startindex, rho=rho_start_value)
```

:1234: **Table**

<table border=2>
<caption align="bottom"> SUMMARY OF m1_rudd_day_20k_high_random_final </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 2.55 </td> <td align="right"> 1.55 </td> <td align="right"> 1.65 </td> <td align="right"> 0.099 </td> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 2 </td> <td align="right"> 2 </td> <td align="right"> 0.04 </td> <td align="right"> 0.962 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 5 </td> <td align="right"> 5 </td> <td align="right"> 1.83 </td> <td align="right"> 0.104 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 2 </td> <td align="right"> 2 </td> <td align="right"> 0.13 </td> <td align="right"> 0.877 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 7 </td> <td align="right"> 7 </td> <td align="right"> 0.15 </td> <td align="right"> 0.994 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 587.34 </td> <td align="right"> 1287 </td> <td align="right"> 1.64 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 226.29 </td> <td align="right"> 1415 </td> <td align="right"> 0.29 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 77.39 </td> <td align="right"> 1114 </td> <td align="right"> 0.1 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 1340.94 </td> <td align="right"> 1712 </td> <td align="right"> 10.31 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.78 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 201716 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 78% </td> </tr>
   <a name=tab.gam></a>
</table>

:chart: **Plots**

![Seiche_models](/outputs/plots/m1_rudd_day_20k_high_random_final.png "Seiche_models")


## <a name="headSpeciesRudd"></a>![Image name](/outputs/icons/rudd_body_1.png) [:waning_crescent_moon:](#head1)

This dataset fits better with a log normal distribution. In general, any of the models fitted showed both poor residuals patterns and edf/k´ ratios. The T-distribution did not result in a better model fit<br />
The k20 model is preferred to models with higher k-values, mixed (lower or higher for the random terms) or not according to Chi-square test and AIC values.<br />
Although the k100 model has higher <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a> than the k20 model, in no case did increasing the k-value improve the edf/k` ratio or the k-index. Along with this, the k100 model could not be re-fitted with fixed degrees of freedom, the significance of the main-effects terms were consistent with the k20 model as well as the general patterns besides the level of wiggliness.<br />

``` r
AIC_table <- AIC(m1_rudd_night_20k_high_global,m1_rudd_night_20k_high_random,m1_rudd_night_20k,m1_rudd_night_10k)%>%
             rownames_to_column(var= "Model")%>%
             mutate(data_source = rep(c("data_rudd_night"), each=1))%>%
             group_by(data_source)%>%
             mutate(deltaAIC = AIC - min(AIC))%>%
             ungroup()%>%
             dplyr::select(-data_source)%>%
             mutate_at(.vars = vars(df,AIC, deltaAIC),
                       .funs = funs(round,.args = list(digits=0)))

kable(AIC_table)
```
|Model                         |  df|   AIC| deltaAIC|
|:-----------------------------|---:|-----:|--------:|
|m1_rudd_night_20k_high_global | 406| 30571|     2119|
|m1_rudd_night_20k_high_random | 682| 28849|      398|
|m1_rudd_night_20k             | 719| 28451|        0|
|m1_rudd_night_10k             | 411| 30087|     1636|


:balance_scale: **Model**

``` r

:balance_scale: **Model**

``` r
m1s_rudd_night_20k <- bam(det_depth ~ s(seasonal_depth, k=20, bs="cr") +
                                      s(amplitude, k=20, bs="cr") +
                                      s(mean_gradient, k=20, bs="cr") +
                                      s(seasonal_depth, fishid, bs="fs", k=20, m=1) +
                                      s(amplitude, fishid, bs="fs", k=20, m=1) +
                                      s(mean_gradient, fishid, bs="fs", k=20, m=1) +
                                      s(dets_ts, k=20, bs="cr") +
                                      s(fishid, dets_ts, bs="fs", k=20, m=1),
                                      data=data_rudd_night,
                                      family=gaussian(link="log"), discrete=TRUE,
                                      nthreads=10, cluster=10, gc.level=0)

rho_start_value <- start_value_rho(m1s_rudd_night_20k, plot=TRUE)

m1_rudd_night_20k <- bam(det_depth ~ s(seasonal_depth, k=20, bs="cr") +
                                     s(amplitude, k=20, bs="cr") +
                                     s(mean_gradient, k=20, bs="cr") +
                                     s(seasonal_depth, fishid, bs="fs", k=20, m=1) +
                                     s(amplitude, fishid, bs="fs", k=20, m=1) +
                                     s(mean_gradient, fishid, bs="fs", k=20, m=1) +
                                     s(dets_ts, k=20, bs="cr") +
                                     s(fishid, dets_ts, bs="fs", k=20, m=1),
                                     data=data_rudd_night,
                                     family=gaussian(link="log"), discrete=TRUE,
                                     nthreads=10, cluster=10, gc.level=0,
                                     AR.start=startindex, rho=rho_start_value)
```
:1234: **Table**

<table border=0>
<caption align="bottom"> SUMMARY OF m1_rudd_night_20k </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 0.67 </td> <td align="right"> 0.51 </td> <td align="right"> 1.31 </td> <td align="right"> 0.189 </td> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 15.71 </td> <td align="right"> 16.47 </td> <td align="right"> 1.67 </td> <th align="right"> 0.026 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 16.04 </td> <td align="right"> 16.91 </td> <td align="right"> 4.43 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 9.3 </td> <td align="right"> 10.38 </td> <td align="right"> 3.03 </td> <th align="right"> 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 17.79 </td> <td align="right"> 18.17 </td> <td align="right"> 1.28 </td> <td align="right"> 0.28 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 139.04 </td> <td align="right"> 245 </td> <td align="right"> 7.38 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 184.06 </td> <td align="right"> 292 </td> <td align="right"> 5.76 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 162.64 </td> <td align="right"> 298 </td> <td align="right"> 2.46 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 165.86 </td> <td align="right"> 255 </td> <td align="right"> 7.44 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.45 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 28451 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 45.6% </td> </tr>
   <a name=tab.gam></a>
</table>

:chart: **Plots**

![Seiche_models](/outputs/plots/m1_rudd_night_20k.png "Seiche_models")


Re-fit the model with fixed degreees of freedom

:balance_scale: **Final model**
``` r
rho_start_value <- start_value_rho(m1s_rudd_night_20k, plot=TRUE)

f.df <- round(summary(m1_rudd_night_20k)$edf)+1 # get edf per smooth
f.df <- pmax(f.df,3)                            # minimum basis dimension is 3

m1_rudd_night_20k_final <- bam(det_depth ~ s(seasonal_depth, k=f.df[1], fx=TRUE) +
                                           s(amplitude, k=f.df[2], fx=TRUE) +
                                           s(mean_gradient, k=f.df[3], fx=TRUE) +
                                           s(seasonal_depth, fishid, bs="fs", k=f.df[4], m=1) +
                                           s(amplitude, fishid, bs="fs", k=f.df[5], m=1) +
                                           s(mean_gradient, fishid, bs="fs", k=f.df[6], m=1) +
                                           s(dets_ts, k=f.df[7], fx=TRUE) +
                                           s(fishid, dets_ts, bs="fs", k=f.df[8], m=1),
                                           data=data_rudd_night,
                                           family=gaussian(link="log"), discrete=TRUE,
                                           nthreads=10, cluster=10, gc.level=0,
                                           AR.start=startindex, rho=rho_start_value)
```

:1234: **Table**

<table border=2>
<caption align="bottom"> SUMMARY OF m1_rudd_night_20k_final </caption>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Parametric&space;\&space;coefficients}" title="\mathbf{Parametric \ coefficients}" /></a> </td> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Estimate}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Estimate}" title="\mathbf{Estimate}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Std.&space;Error}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Std.&space;Error}" title="\mathbf{Std. Error}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{t}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{t}" title="\mathbf{t}" /></a> </th> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a> </td> <td align="right"> 0.88 </td> <td align="right"> 1.01 </td> <td align="right"> 0.87 </td> <td align="right"> 0.39 </td> </tr>
   <tr> <td align="left"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{edf}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{edf}" title="\boldsymbol{edf}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{Ref.df}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{Ref.df}" title="\boldsymbol{Ref.df}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{F}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{F}" title="\mathbf{F}" /></a> </td> <td align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p\,&space;\mathbf{-value}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{p\,&space;\mathbf{-value}}" title="\boldsymbol{p\, \mathbf{-value}}" /></a> </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a> </td> <td align="right"> 16 </td> <td align="right"> 16 </td> <td align="right"> 1.36 </td> <td align="right"> 0.15 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a> </td> <td align="right"> 16 </td> <td align="right"> 16 </td> <td align="right"> 2.52 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a> </td> <td align="right"> 9 </td> <td align="right"> 9 </td> <td align="right"> 1.07 </td> <td align="right"> 0.38 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a> </td> <td align="right"> 18 </td> <td align="right"> 18 </td> <td align="right"> 0.3 </td> <td align="right"> 0.998 </td> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 636.2 </td> <td align="right"> 1233 </td> <td align="right"> 2.67 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 856.2 </td> <td align="right"> 2170 </td> <td align="right"> 1.32 </td> <th align="right"> &lt; 0.001 </th> </tr>
 <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 347.1 </td> <td align="right"> 2178 </td> <td align="right"> 0.37 </th> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <td> <a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a> </td> <td align="right"> 1300.8 </td> <td align="right"> 1562 </td> <td align="right"> 10.48 </td> <th align="right"> &lt; 0.001 </th> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a> </th> <td align="right"> 0.63 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a> </th> <td align="right"> 1710 </td> </tr>
  <tr> <th align="right"> <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a> </td> <td align="right"> 62.8% </td> </tr>
   <a name=tab.gam></a>
</table>

:chart: **Plots**

![Seiche_models](/outputs/plots/m1_rudd_night_20k_final.png "Seiche_models")


# <a name="headafinaltable"></a>Final models table [:page_facing_up:](#headindex)



|    <br>                            |           <br><a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Pike}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Pike}" title="\mathbf{Pike}" /></a>           |                                |       <br><a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Wels\&space;catfish}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Wels\&space;catfish}" title="\mathbf{Wels\ catfish}" /></a>       |                                |           <br><a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Tench}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Tench}" title="\mathbf{Tench}" /></a>          |                                |           <br><a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Rudd}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Rudd}" title="\mathbf{Rudd}" /></a>           |                                |
|------------------------------------|:----------------------------:|:------------------------------:|:----------------------------:|:------------------------------:|:----------------------------:|:------------------------------:|:----------------------------:|:------------------------------:|
|    <br>                            |    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{Day}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathrm{Day}" title="\mathrm{Day}" /></a><br>            |    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{Night}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathrm{Night}" title="\mathrm{Night}" /></a><br>            |    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{Day}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathrm{Day}" title="\mathrm{Day}" /></a><br>            |    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{Night}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathrm{Night}" title="\mathrm{Night}" /></a><br>            |    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{Day}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathrm{Day}" title="\mathrm{Day}" /></a><br>            |    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{Night}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathrm{Night}" title="\mathrm{Night}" /></a><br>            |    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{Day}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathrm{Day}" title="\mathrm{Day}" /></a><br>            |    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{Night}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathrm{Night}" title="\mathrm{Night}" /></a><br>            |
|    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;\mathbf{Parametric&space;\&space;coefficients}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\large&space;\mathbf{Parametric&space;\&space;coefficients}" title="\large \mathbf{Parametric \ coefficients}" /></a></a>  |                              |                                |    <br>                      |    <br>                        |    <br>                      |    <br>                        |    <br>                      |    <br>                        |
|    <br><a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;(&space;Intercept&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\left&space;(&space;Intercept&space;\right&space;)" title="\left ( Intercept \right )" /></a>                 |           <br>5.64           |            <br>4.36            |           <b> <p>1.59**      |           <b> <p>1.64**        |           <br>3.55           |            <br>0.92            |           <br>2.55           |            <br>0.88            |
|    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Smooth&space;\&space;terms}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{Smooth&space;\&space;terms}" title="\mathbf{Smooth \ terms}" /></a> </td>             |                              |                                |    <br>                      |    <br>                        |    <br>                      |    <br>                        |    <br>                      |    <br>                        |
|    <br><a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{seasonal&space;\&space;depth}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{seasonal \ depth} \right )" /></a>           |             <br>2            |              <br>8             |             <br>2            |           <b> <p>17*           |            <br>6             |              <br>8             |             <br>2            |             <br>16             |
|    <br><a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{amplitude}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{amplitude} \right )" /></a>                |            <br>2         |             <b> <p>8***        |             <br>2            |             <br>13             |             <br>5            |              <b> <p>5**        |             <br>5            |             <b> <p>16***       |
|    <br><a href="https://www.codecogs.com/eqnedit.php?latex=\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathit{s}\left&space;(&space;\mathrm{mean&space;\&space;gradient}&space;\right&space;)" title="\mathit{s}\left ( \mathrm{mean \ gradient} \right )" /></a>            |             <br>2            |              <br>3             |           <b> <p>4***        |              <br>2             |             <b> <p>6*        |            <b> <p>6***         |             <br>2            |              <br>9             |
|    <br><a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)" title="s\left( \mathrm{time} \right )" /></a>                     |             <br>2            |              <br>8             |             <br>2            |             <br>17             |            <br>9             |              <br>9             |             <br>7            |             <br>18             |
|    <br><a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{seasonal\&space;depth}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{seasonal\ depth} \right )_{Fish\ ID}" /></a>    |            <b> <p>904.9***   |            <b> <p>709***       |          <b> <p>572.3***     |     <b> <p>736.2***            |    <b> <p>766.54***          |     <b> <p>629.4***            |    <b> <p>587.34***          |     <b> <p>636.2***            |
|    <br><a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{amplitude}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{amplitude} \right )_{Fish\ ID}" /></a>         |            <b> <p>378.1***   |           <b> <p>454.7***      |          <b> <p>270.68***    |           <b> <p>297.5***      |         <b> <p>110.64***     |            <b> <p>448***       |          <b> <p>226.29***    |           <b> <p>856.2***      |
|    <br><a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{mean\&space;gradient}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{mean\ gradient} \right )_{Fish\ ID}" /></a>     |            <b> <p>139.1***   |           <b> <p>276.7***      |          <b> <p>63.25***     |           <b> <p>209.4***      |         <b> <p>84.55***      |           <b> <p>129.5***      |          <b> <p>77.39***     |           <b> <p>347.1***      |
|    <br><a href="https://www.codecogs.com/eqnedit.php?latex=s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?s\left(&space;\mathrm{time}&space;\right&space;)_{Fish\&space;ID}" title="s\left( \mathrm{time} \right )_{Fish\ ID}" /></a>              |           <b> <p>1787***   |           <b> <p>1012.6***     |       <b> <p>2871.52***      |          <b> <p>1369.6***      |           <b> <p>1464.94***  |           <b> <p>1139.9***     |         <b> <p>1340.94***    |           <b> <p>1300.8***     |
|    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{N}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{N}" title="\boldsymbol{N}" /></a>                        |           <br>12             |            <br>12              |           <br>15             |            <br>15              |           <br>19             |             <br>19             |           <br>15             |            <br>15              |
|    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Adjusted\&space;\boldsymbol{R^{2}}}" title="\mathbf{Adjusted\ \boldsymbol{R^{2}}}" /></a>              |           <br>0.92           |            <br>0.85            |           <br>0.78           |            <br>0.76            |           <br>0.63           |             <br>0.6            |           <br>0.78           |            <br>0.63           |
|    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{AIC}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{AIC}" title="\boldsymbol{AIC}" /></a>                      |           <br>204,822           |           <br>53,694           |          <br>659,663         |           <br>330,326          |          <br>864,749        |           <br>188,051          |          <br>201,716         |            <br>1710            |
|    <b> <p><a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{Deviance\&space;explained}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{Deviance\&space;explained}" title="\mathbf{Deviance\ explained}" /></a>       |          <br>91.7%          |           <br>75.40%           |          <br>78.40%          |           <br>68.10%           |          <br>55.60%          |           <br>59.30%           |            <br>78%           |           <br>62.80%           |

<a href="https://www.codecogs.com/eqnedit.php?latex=\small&space;\mathrm{*&space;\mathit{p}<0.05,&space;**&space;\mathit{p}<0.01,&space;***&space;\mathit{p}<0.001}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\small&space;\mathrm{*&space;\mathit{p}<0.05,&space;**&space;\mathit{p}<0.01,&space;***&space;\mathit{p}<0.001}" title="\small \mathrm{* \mathit{p}<0.05, ** \mathit{p}<0.01, *** \mathit{p}<0.001}" /></a>

# <a name="headaltmodels"></a> Alternative model structures [:page_facing_up:](#headindex)

### <a name="head6"></a>:balance_scale: _Model 2: Group-level smoothers with same wiggliness (Random factor-smoother interaction)_

Depth is modelled as a function of four smoother terms:
- 0 *main-effects smoothers* (not valid if we want to estimate average sesonal depth, mean gradient or amplitude).
- 0 group-specific *random-effect* smoothers (i.e., without fish-specific intercepts).
- 4 group-level *random factor-smoother interactions* (i.e., for all _Fish ID_ with `bs="fs"`):
  - A penalty to shrink these smoothers toward zero.
  - A common smoothing parameter for all _Fish ID_ smoothers with same wiggliness (<a href="https://www.codecogs.com/eqnedit.php?latex=\simeq" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\simeq" title="\simeq" /></a> functional responses for all _Fish ID_)
  - Individual shapes of the smooth terms not related.
``` r
m2 <- bam(det_depth ~ s(seasonal_depth, fishid, bs="fs", k=10, m=1) +
                      s(amplitude, fishid, bs="fs", k=10, m=1) +
                      s(mean_gradient, fishid, bs="fs", k=10, m=1) +
                      s(dets_ts, fishid, bs="fs", k=10, m=1),
                      data=data_wels_day, method="REML", family="gaussian",
                      discrete=TRUE, nthreads=10, cluster=10, gc.level=0,
                      AR.start=startindex, rho=rho_start_value)
```

### <a name="head6"></a>:balance_scale: _Model 3: Main-effects smoothers + Group-level smoothers with different wiggliness_

Now, using the argument `bs="re"`, the random effects are not being smoothed but just _"exploiting the link between smooths and random effects to treat random effects as smooths"_ (see `?mgcv::smooth.construct.re.smooth.spec`). Thus, for each level factor in _fishid_, an independent normal prior with common variance is assumed.<br />
Depth is modelled as a function of nine smoother terms:
- 4 *main-effects* smoothers (with cubic splines) for all observations of each covariate (individual effects with individual penalty terms).
- 1 group-specific *random-effects* smoother (i.e., fish-specific intercepts with `bs="re"`).
- 4 group-level *factor-smoothers* (i.e., _Fish ID_-specific with `by=fishid`):
  - A penalty to shrink these smoothers toward zero.
  - A different smoothing parameter for each _Fish ID_ smoother, each with different wiggliness (<a href="https://www.codecogs.com/eqnedit.php?latex=\neq" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\neq" title="\neq" /></a> functional responses for each _Fish ID_)
  - Different shapes of the smooth terms (~inter-individual variation in responses, t2).
``` r
m3 <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
		                  s(amplitude, bs="cr", k=10) +
                      s(mean_gradient, bs="cr", k=10) +
                      s(dets_ts, bs="cr", k=10) +
                      s(seasonal_depth, by=fishid, bs="cr", k=10, m=1) +
                      s(amplitude, by=fishid, bs="cr", k=10, m=1) +
                      s(mean_gradient, by=fishid, bs="cr", k=10, m=1) +
                      s(dets_ts, by=fishid, bs="cr", k=10, m=1) +
                      s(fishid, bs="re", k=14),
                      data=data_wels_day, method="REML", family="gaussian",
                      discrete=TRUE, nthreads=10, cluster=10, gc.level=0,
                      AR.start=startindex, rho=rho_start_value)
```

### <a name="head6"></a>:balance_scale: _Model 4: Group-level smoothers with different wiggliness_

Depth is modelled as a function of five smoother terms:
- 0 *main-effects smoothers* (not valid if we want to estimate average sesonal depth, mean gradient or amplitude).
- 1 group-specific *random-effects* smoother (i.e., fish-specific intercepts with `bs="re"`).
- 4 group-level *factor-smoothers* (i.e., _Fish ID_-specific with `by=fishid`):
  - A penalty to shrink these smoothers toward zero.
  - A different smoothing parameter for each _Fish ID_ smoother, each with different wiggliness (<a href="https://www.codecogs.com/eqnedit.php?latex=\neq" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\neq" title="\neq" /></a> functional responses for each _Fish ID_)
  - Different shapes of the smooth terms (~inter-individual variation in responses, t2).
``` r
m4 <- bam(det_depth ~ s(seasonal_depth, by=fishid, bs="cr", k=10, m=1) +
                      s(amplitude, by=fishid, bs="cr", k=10, m=1) +
                      s(mean_gradient, by=fishid, bs="cr", k=10, m=1) +
                      s(dets_ts, by=fishid, bs="cr", k=10, m=1) +
                      s(fishid, bs="re", k=14),
                      data=data_wels_day, method="REML", family="gaussian",
                      discrete=TRUE, nthreads=10, cluster=10, gc.level=0,
                      AR.start=startindex, rho=rho_start_value)
```
:information_source: _Note that _s(fishid, bs="re")_ is included in models 3 and 4 for an intercept as the (_X_ `by` _fishid_) smoother does not include one._ In both cases, it would make more sense to include the grouping variable _species_ instead of _fishid_.

### <a name="head6"></a>:balance_scale: _Model 5: Main-effects smoothers + Random-effects smoothers (Intercept)_

Depth is modelled as a function of five smoother terms:
- 4 *main-effects* smoothers (cubic splines) for all observations of each covariate.
- 1 group-specific *random-effects* smoother (i.e., fish-specific intercepts with `bs="re"`).
``` r
m5 <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
		                  s(amplitude, bs="cr", k=10) +
                      s(mean_gradient, bs="cr", k=10) +
                      s(dets_ts, bs="cr", k=10) +
                      s(fishid, bs="re", k=14),
                      data=data_wels_day, method="REML", family="gaussian",
                      discrete=TRUE, nthreads=10, cluster=10, gc.level=0,
                      AR.start=startindex, rho=rho_start_value)
```

### <a name="head6"></a>:balance_scale: _Model 6: Main-effects smoothers + Random-effects smoothers (Intercept + slope)_

Depth is modelled as a function of six smoother terms:
- 4 *main-effects* smoothers (cubic splines) for all observations of each covariate.
- 1 group-specific *random-effects* smoother (i.e., fish-specific intercepts with `bs="re"`).
- 3 group-specific *random-effects* smoothers (i.e., fish-specific slopes with `bs="re"`).
``` r
m6 <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
		                  s(amplitude, bs="cr", k=10) +
                      s(mean_gradient, bs="cr", k=10) +
                      s(dets_ts, bs="cr", k=10) +
                      s(seasonal_depth, fishid, bs="re", k=14) +
                      s(amplitude, fishid, bs="re", k=14) +
                      s(mean_gradient, fishid, bs="re", k=14) +
                      s(fishid, dets_ts, bs="re", k=14) +
                      s(fishid, bs="re", k=14),
                      data=data_wels_day, method="REML", family="gaussian",
                      discrete=TRUE, nthreads=10, cluster=10, gc.level=0,
                      AR.start=startindex, rho=rho_start_value)
```

### <a name="head6"></a>:balance_scale: _Models with sesonal dynamics_

Consider any of the above models for including a seasonal components. For example, for the selected **Model 1**:
- Include cyclic cubic splines (`bs="cc"`) for the daily dynamics (day-to-day effects cyclic variation in depth).
- Specify start and end points of cycles (1 day = 1440 minutes).

Depth is modelled as a function of nine smoother terms:
- 4 *main-effects* smoothers for all observations of each covariate (individual effects have a common penalty term):
  - With cubic splines (3).
  - With  cyclic cubic spline (1).
- 1 group-specific *random-effects* smoother (i.e., fish-specific slopes with `bs="re"`) for daily dynamics.
- 4 group-level *random factor-smoother interactions* (i.e., for all _Fish ID_ with `bs="fs"`):
  - A penalty to shrink these smoothers toward zero.
  - A common smoothing parameter for all _Fish ID_ smoothers with same wiggliness (<a href="https://www.codecogs.com/eqnedit.php?latex=\simeq" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\simeq" title="\simeq" /></a> functional responses for all _Fish ID_)
  - Different shapes of the smooth terms (~inter-individual variation in responses, t2).
  - One incorporating sesonal dynamics:
    - With the `knots` argument for start and end cycles.
    - With cyclic cubic spline (1).
``` r
m1_seasonal <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
                               s(amplitude, bs="cr", k=10) +
                               s(mean_gradient, bs="cr", k=10) +
                               s(dets_ts, k=10, bs='cc') +
                               s(seasonal_depth, fishid, bs="fs", k=10, m=1) +
                               s(amplitude, fishid, bs="fs", k=10, m=1) +
                               s(mean_gradient, fishid, bs="fs", k=10, m=1) +
                               s(fishid, day, bs="re") +
                               s(dets_ts, fishid, k=10, bs="fs", xt=list(bs="cc")), knots=list(day=c(0,1440)),
                               data=data_wels_day, method="REML", family="gaussian",
                               discrete=TRUE, nthreads=10, cluster=10, gc.level=0,
                               AR.start=startindex, rho=rho_start_value)
```

Based on model 4:
``` r
m4_seasonal <- bam(det_depth ~ s(seasonal_depth, by=fishid, k=10, bs="cc") +
                               s(amplitude, by=fishid, k=10, bs="cc") +
                               s(mean_gradient, by=fishid, k=10, bs="cc") +
                               s(dets_ts, by=fishid, k=10, bs="cc") +
                               s(fishid, bs="re") +
                               s(fishid, day, bs="re") + knots=list(day=c(0,1440)),
                               data=data_wels_day, method="REML", family="gaussian",
                               discrete=TRUE, nthreads=10, cluster=10, gc.level=0,
                               AR.start=startindex, rho=rho_start_value)
```

# <a name="headissuesandimprovements"></a> Issues and further improvements [:page_facing_up:](#headindex)

For each subset of data, shifting in k-value from 10 to 20 to 50 to 100 increased <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a> but as the number of basis functions increased, models were less likely to be re-fitted with fixed degrees of freedom at k=50-100. Therefore, for each specific data set we kept the model with the upper limit of flexibility allowed provided a balance between a reasonable k-value, AIC, <a href="https://www.codecogs.com/eqnedit.php?latex=R&space;^{2}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?R&space;^{2}" title="R ^{2}" /></a> and Chi-square test results. In most cases, random smooths with ten basis functions (k=10) were not sufficient for modelling some of the data, otherwise setting k = 20 provided more reliable estimates.<br />
In general, the results of modelling using different basis functions evidence that the inclusion of individual random smooth terms for each covariate makes the model appear to never get an optimal fit as more basis dimensions are needed for the higher order of the random terms. In other words, increasing k is fine for estimating the main smooth effects but for the random smooth effects, the higher variability would be better modelled with regularly sparsed chunks of data. This, in part, is atributible to the fact that without a seasonal component complexity of the already complex time-series trend increases thus complicating it to be approximated by a relatively small number of basis functions.<br />
One solution would be increasing the k-value to theoretically higher k-values (for all or only the random terms) but this confronts two problems, either the model can`t be fitted with so many basis functions or, if fitted, it´ll never be re-fitted with un-penalized terms. In both cases, what we see are patterns of excessive wiggliness in the random terms. This also tends to change the estimation of the main-effects terms such as p-values can shift from significance to non-significance and viceversa.<br />
I highlight this and some other issues with possible solutions that could be further worked out in a future revision.

### <a name="head6"></a>:wrench: _Lack of continutiy in the data_

- Each subset consists of discontinuous data separated into day-night. For a specific day, this causes an individual`s measurements to go from day to day without the intermediate night values.
- A way of solving this would be including the factor variable _diel period_ in the models (for day-night differences), for example using tensor product interactions (`te()`) and including _diel_period_ using the argument `"by"`.
- Discontinuity is also due to missing data points within diel periods. This could be solved by computing the missing cases before fitting the model.
- Examples of model formulas including _diel_period_:

``` r
m4_seasonal_diel <- bam(det_depth ~ diel_period +
                                    te(day, seasonal_depth, by=diel_period, bs=c("cc","tp"), k=c(10,10), m=1) +
                                    te(day, amplitude, by=diel_period, bs=c("cc","tp"), k=c(10,10), m=1) +
                                    te(day, mean_gradient, by=diel_period, bs=c("cc","tp"), k=c(10,10), m=1) +
                                    te(day, dets_ts, by=diel_period, bs=c("cc","tp"), k=c(10,10), m=1) +
                                    s(fishid, bs="re") +
                                    s(fishid, day, bs="re") + knots=list(day=c(0,1440)),
                                    data=data_wels, method="REML", family="gaussian",
                                    discrete=TRUE, nthreads=10, cluster=10, gc.level=0,
                                    AR.start=startindex, rho=rho_start_value)
```
:information_source: _Note that I use here `"tp"` splines for the main-effects smoothers instead of `"cr"` and cyclic cubic splines `"cc"` for the random smoothers._


### <a name="head6"></a>:wrench: _Lack of sesonality_

- Even though accounting for sesonality is out of the scope of this study, including a seasonal component such as shown in **Models with sesonal dynamics** above would help improve the fitting and predictability of the model.
- Could this be fixed using adaptive splines (`bs="ad"`)?
  - Uses a weighted penalty matrix with the weights varying smoothly with the covariate to accommodate the varying wiggliness.
  - Uses more smoothness parameters (`m=5` by default) but assumes the same degree of smoothness over the range of the covariate.

### <a name="head6"></a>:wrench: _Including logger sites for predicting changes in thermocline?_

- Similarly as above, including the two logger sites with the argument `"by"` can be used to model variability between each section of the lake (eastern-western).
- If we model depth as a function of time and logger site we may include an interaction using a smooth _tensor product interaction_ between logger site and time. Modifying **Model 3** as follows:

``` r
m3_logger_site <- bam(det_depth ~ s(seasonal_depth, bs="cr", k=10) +
		                              s(amplitude, bs="cr", k=10) +
                                  s(mean_gradient, bs="cr", k=10) +
                                  s(dets_ts, bs="cr", k=10) +
                                  s(seasonal_depth, by=logger_site, bs="cr", k=10, m=1) +
                                  s(amplitude, by=logger_site, bs="cr", k=10, m=1) +
                                  s(mean_gradient, by=logger_site, bs="cr", k=10, m=1) +
                                  s(dets_ts, by=logger_site, bs="cr", k=10, m=1) +
                                  s(fishid, bs="re", k=14),
                                  data=data_wels_day, method="REML", family="gaussian",
                                  discrete=TRUE, nthreads=10, cluster=10, gc.level=0,
                                  AR.start=startindex, rho=rho_start_value)
```

### <a name="head6"></a>:wrench: _Variable selection and parameters estimation of the main and random smooth effects_

- All models are fitted to the same global formula (**Model 1**) so there`s not possibility for variable selection. This can be an issue when non-significant main-effects smoothers tend to bias estimates and p-values of other smooth terms.
- Variable selection can be done using more reliable p-values estimation of the main smooth effects by using the argument `select=TRUE` that penalizes all model terms to zero effect. Even if we don`t use this argument, the fitted models can have potentially higher deviance from including non-significant terms.
- Rather using `bs="fs"` we could use the multivariate equivalent of the _factor-smoother interaction_ `t2(full=TRUE)` that appears to estimate the _main-effects smoothers_ more accurately when _group-level smoothers_ are present.
- The problem with the former two approaches is that the computational cost is very high.


# <a name="headreferences"></a> References [:page_facing_up:](#headindex)

- **Marra, G., & Wood, S. N. (2011)**. Practical variable selection for generalized additive models. Computational Statistics & Data Analysis. 2011;55(7):2372–2387. https://doi.org/10.1016/j.csda.2011.02.004
- **Pedersen, E. J., Miller, D. L., Simpson, G. L., & Ross, N. (2019)**. Hierarchical generalized additive models in ecology: an introduction with mgcv. PeerJ, 7, e6876. https://doi.org/10.7717/peerj.6876
- **van Rij, J., Wieling, M., Baayen, R., van Rijn, H. (2017)**. itsadug: Interpreting Time Series and Autocorrelated Data Using GAMMs. R package version 2.3
- **Wood, S. N. (2001)**. mgcv: GAMs and Generalized Ridge Regression for R. R News 1(2):20-25
- **Wood, S. N. (2017)**. Generalized additive models: an introduction with R. Second Edition. Boco Raton: CRC Press.
- **Wood, S. N., Goude, Y., & Shaw, S. (2015)**. Generalized additive models for large datasets. Journal of the Royal Statistical Society, Series C 64(1): 139-155. http://dx.doi.org/10.1111/rssc.12068
- **Wood, S. N., Pya, N., & Saefken, B. (2016)**. Smoothing parameter and model selection for general smooth models. Journal of the American Statistical Association, 111:516, 1548-1563. https://doi.org/10.1080/01621459.2016.1180986


