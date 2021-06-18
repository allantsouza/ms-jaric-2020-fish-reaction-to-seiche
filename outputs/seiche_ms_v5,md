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
