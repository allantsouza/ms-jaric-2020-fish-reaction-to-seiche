#Libraries
library(lme4)
library(lmerTest)
library(MuMIn)
library(sjPlot)
library(bbmle)
library(hier.part)
library(parallel)
library(tictoc)#to check time elapsed
library(optimx)#optmization for lmer
library(xlsx)

#Function to run lmer in multiple cores
f_lmer_mc <- function(data, calls, mc.cores) {
  require(parallel)
  if (is.data.frame(data)) 
    data = replicate(length(calls), data, simplify = F)
  for (i in 1:length(data)) attr(data[[i]], "cll") = calls[i]
  m.list = mclapply(data, function(i) eval(parse(text = attr(i, "cll"))), 
                    mc.cores = mc.cores)
  return(m.list)
}



# Loading the data ####
fish_raw <- read_csv(file = "data/raw/fishIDs.csv", col_types = "ccdc") %>%
  mutate(data_path = here("data/products/fish",paste0(tag_sn, ".csv")))

detections <- fish_raw %>%
  filter(file.exists(data_path)) %>%
  pull(data_path) %>%
  # read in all the files, appending the path before the filename
  map(~ read_csv(.)) %>% 
  reduce(rbind) %>%
  inner_join(fish_raw[,c("tag_sn","fishid", "species")])


#Models 2020####
#Pike ####
#Models run - Pike ####
#random intercept and random slope per fish id - Day (Seiche only)
tic('It took this amount of time to run the model')
mdl_random_slo_int_day_pike <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + 
                                 det_therm_strength + lake_therm_depth_smoothed_center + 
                                 det_therm_deviation_center + 
                                 (1 + lake_therm_thickness_smoothed + 
                                    det_therm_strength + 
                                    lake_therm_depth_smoothed_center +
                                    det_therm_deviation_center | fishid),
                               data = detections %>% 
                                 filter(species == "pike" & 
                                          diel_period == 'day' & 
                                          is_valid_seiche == TRUE), 
                               REML = T, 
                               lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()

#random intercept and random slope per fish id - Night (Seiche only)
tic('It took this amount of time to run the model')
mdl_random_slo_int_night_pike <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + 
                                 det_therm_strength + lake_therm_depth_smoothed_center + 
                                 det_therm_deviation_center + 
                                 (1 + lake_therm_thickness_smoothed + 
                                    det_therm_strength + 
                                    lake_therm_depth_smoothed_center + 
                                    det_therm_deviation_center | fishid),
                               data = detections %>% 
                                 filter(species == "pike" & 
                                          diel_period == 'night' & 
                                          is_valid_seiche == TRUE), 
                               REML = T, 
                               lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()

#Summary of the model results - Pike ####
tab_model(mdl_random_slo_int_day_pike, mdl_random_slo_int_night_pike, dv.labels = c("Day","Night"), 
          title='Pike (Seiche only)',
          pred.labels = c('Intercept', 'Thermocline thickness (m)',
                          'Thermocline strength (ºC)',
                          'Thermocline depth (m)',
                          'Seiche strength'))#sigma² = withing group variance; tau00 = between group variance; ICC = proportion of variance explained for each group level

#Fixed effects summary - Pike - Day ####
params_mld_pike_day <- parameters::model_parameters(mdl_random_slo_int_day_pike)

mdl_pike_day_fixed <- tibble(species = "pike", 
                             diel_period = "day",
                             parameter = c(params_mld_pike_day$Parameter[1],
                                           params_mld_pike_day$Parameter[2],
                                           params_mld_pike_day$Parameter[3],
                                           params_mld_pike_day$Parameter[4],
                                           params_mld_pike_day$Parameter[5]),
                             coefficient = c(params_mld_pike_day$Coefficient[1],
                                             params_mld_pike_day$Coefficient[2],
                                             params_mld_pike_day$Coefficient[3],
                                             params_mld_pike_day$Coefficient[4],
                                             params_mld_pike_day$Coefficient[5]),
                             se = c(params_mld_pike_day$SE[1],
                                    params_mld_pike_day$SE[2],
                                    params_mld_pike_day$SE[3],
                                    params_mld_pike_day$SE[4],
                                    params_mld_pike_day$SE[5]),
                             ci_05 = c(params_mld_pike_day$CI_low[1],
                                       params_mld_pike_day$CI_low[2],
                                       params_mld_pike_day$CI_low[3],
                                       params_mld_pike_day$CI_low[4],
                                       params_mld_pike_day$CI_low[5]),
                             ci_95 = c(params_mld_pike_day$CI_high[1],
                                       params_mld_pike_day$CI_high[2],
                                       params_mld_pike_day$CI_high[3],
                                       params_mld_pike_day$CI_high[4],
                                       params_mld_pike_day$CI_high[5]),
                             t_value = c(params_mld_pike_day$t[1],
                                         params_mld_pike_day$t[2],
                                         params_mld_pike_day$t[3],
                                         params_mld_pike_day$t[4],
                                         params_mld_pike_day$t[5]),
                             p_value = c(params_mld_pike_day$p[1],
                                         params_mld_pike_day$p[2],
                                         params_mld_pike_day$p[3],
                                         params_mld_pike_day$p[4],
                                         params_mld_pike_day$p[5]))

#Random effects summary - Pike - Day ####
rdm_mdl_pike_day <- parameters::random_parameters(mdl_random_slo_int_day_pike)

mdl_pike_day_random <- tibble(species = "pike",
                              diel_period = "day",
                              description = c(rdm_mdl_pike_day$Description[1],
                                              rdm_mdl_pike_day$Description[2],
                                              rdm_mdl_pike_day$Description[3],
                                              rdm_mdl_pike_day$Description[4],
                                              rdm_mdl_pike_day$Description[5],
                                              rdm_mdl_pike_day$Description[6],
                                              rdm_mdl_pike_day$Description[7],
                                              rdm_mdl_pike_day$Description[8],
                                              rdm_mdl_pike_day$Description[9],
                                              rdm_mdl_pike_day$Description[10],
                                              rdm_mdl_pike_day$Description[11],
                                              rdm_mdl_pike_day$Description[12],
                                              rdm_mdl_pike_day$Description[13],
                                              rdm_mdl_pike_day$Description[14],
                                              rdm_mdl_pike_day$Description[15],
                                              rdm_mdl_pike_day$Description[16],
                                              rdm_mdl_pike_day$Description[17],
                                              rdm_mdl_pike_day$Description[18]),
                              component = c(rdm_mdl_pike_day$Component[1],
                                            rdm_mdl_pike_day$Component[2],
                                            rdm_mdl_pike_day$Component[3],
                                            rdm_mdl_pike_day$Component[4],
                                            rdm_mdl_pike_day$Component[5],
                                            rdm_mdl_pike_day$Component[6],
                                            rdm_mdl_pike_day$Component[7],
                                            rdm_mdl_pike_day$Component[8],
                                            rdm_mdl_pike_day$Component[9],
                                            rdm_mdl_pike_day$Component[10],
                                            rdm_mdl_pike_day$Component[11],
                                            rdm_mdl_pike_day$Component[12],
                                            rdm_mdl_pike_day$Component[13],
                                            rdm_mdl_pike_day$Component[14],
                                            rdm_mdl_pike_day$Component[15],
                                            rdm_mdl_pike_day$Component[16],
                                            rdm_mdl_pike_day$Component[17],
                                            rdm_mdl_pike_day$Component[18]),
                              type = c(rdm_mdl_pike_day$Type[1],
                                       rdm_mdl_pike_day$Type[2],
                                       rdm_mdl_pike_day$Type[3],
                                       rdm_mdl_pike_day$Type[4],
                                       rdm_mdl_pike_day$Type[5],
                                       rdm_mdl_pike_day$Type[6],
                                       rdm_mdl_pike_day$Type[7],
                                       rdm_mdl_pike_day$Type[8],
                                       rdm_mdl_pike_day$Type[9],
                                       rdm_mdl_pike_day$Type[10],
                                       rdm_mdl_pike_day$Type[11],
                                       rdm_mdl_pike_day$Type[12],
                                       rdm_mdl_pike_day$Type[13],
                                       rdm_mdl_pike_day$Type[14],
                                       rdm_mdl_pike_day$Type[15],
                                       rdm_mdl_pike_day$Type[16],
                                       rdm_mdl_pike_day$Type[17],
                                       rdm_mdl_pike_day$Type[18]),
                              term = c(rdm_mdl_pike_day$Term[1],
                                       rdm_mdl_pike_day$Term[2],
                                       rdm_mdl_pike_day$Term[3],
                                       rdm_mdl_pike_day$Term[4],
                                       rdm_mdl_pike_day$Term[5],
                                       rdm_mdl_pike_day$Term[6],
                                       rdm_mdl_pike_day$Term[7],
                                       rdm_mdl_pike_day$Term[8],
                                       rdm_mdl_pike_day$Term[9],
                                       rdm_mdl_pike_day$Term[10],
                                       rdm_mdl_pike_day$Term[11],
                                       rdm_mdl_pike_day$Term[12],
                                       rdm_mdl_pike_day$Term[13],
                                       rdm_mdl_pike_day$Term[14],
                                       rdm_mdl_pike_day$Term[15],
                                       rdm_mdl_pike_day$Term[16],
                                       rdm_mdl_pike_day$Term[17],
                                       rdm_mdl_pike_day$Term[18]),
                              value = c(rdm_mdl_pike_day$Value[1],
                                        rdm_mdl_pike_day$Value[2],
                                        rdm_mdl_pike_day$Value[3],
                                        rdm_mdl_pike_day$Value[4],
                                        rdm_mdl_pike_day$Value[5],
                                        rdm_mdl_pike_day$Value[6],
                                        rdm_mdl_pike_day$Value[7],
                                        rdm_mdl_pike_day$Value[8],
                                        rdm_mdl_pike_day$Value[9],
                                        rdm_mdl_pike_day$Value[10],
                                        rdm_mdl_pike_day$Value[11],
                                        rdm_mdl_pike_day$Value[12],
                                        rdm_mdl_pike_day$Value[13],
                                        rdm_mdl_pike_day$Value[14],
                                        rdm_mdl_pike_day$Value[15],
                                        rdm_mdl_pike_day$Value[16],
                                        rdm_mdl_pike_day$Value[17],
                                        rdm_mdl_pike_day$Value[18]))

#Percentage explained by each variable - Pike - Day ####
#Sum of squares to give the percentage explained by each variable in the model
#Day
mdl_random_slo_int_pike_day_aov <- anova(mdl_random_slo_int_day_pike)
pike_mdl_day_percentage <- tibble(species = "pike",
                                  diel_period = "day",
                                  variable = c('Thermocline thickness (m)',
                                               'Thermocline strength (ºC)',
                                               'Thermocline depth (m)',
                                               'Seiche strength'),
                                  percentage_explained = c(
                                    round(mdl_random_slo_int_pike_day_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_pike_day_aov$`Sum Sq`),2),
                                    round(mdl_random_slo_int_pike_day_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_pike_day_aov$`Sum Sq`),2),
                                    round(mdl_random_slo_int_pike_day_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_pike_day_aov$`Sum Sq`),2),
                                    round(mdl_random_slo_int_pike_day_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_pike_day_aov$`Sum Sq`),2)
                                  )
)

pike_day_mdl_summary <- list(fixed_effects = mdl_pike_day_fixed,
                             percentage_explained = pike_mdl_day_percentage,
                             random_effects = mdl_pike_day_random)

#Night
mdl_random_slo_int_pike_night_aov <-  anova(mdl_random_slo_int_night_pike)
pike_mdl_night_percentage <- tibble(species = "pike",
                                    diel_period = "night",
                                    variable = c('Thermocline thickness (m)',
                                                 'Thermocline strength (ºC)',
                                                 'Thermocline depth (m)',
                                                 'Seiche strength'),
                                    percentage_explained = c(
                                      round(mdl_random_slo_int_pike_night_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_pike_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_pike_night_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_pike_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_pike_night_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_pike_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_pike_night_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_pike_night_aov$`Sum Sq`),2))
)

#pike (day & night)
pike_mld_percentages <- bind_rows(pike_mdl_day_percentage, pike_mdl_night_percentage)

#Fixed effects summary - Pike - Night ####
params_mld_pike_night <- parameters::model_parameters(mdl_random_slo_int_night_pike)

mdl_pike_night_fixed <- tibble(species = "pike", 
                               diel_period = "night",
                               parameter = c(params_mld_pike_night$Parameter[1],
                                             params_mld_pike_night$Parameter[2],
                                             params_mld_pike_night$Parameter[3],
                                             params_mld_pike_night$Parameter[4],
                                             params_mld_pike_night$Parameter[5]),
                               coefficient = c(params_mld_pike_night$Coefficient[1],
                                               params_mld_pike_night$Coefficient[2],
                                               params_mld_pike_night$Coefficient[3],
                                               params_mld_pike_night$Coefficient[4],
                                               params_mld_pike_night$Coefficient[5]),
                               se = c(params_mld_pike_night$SE[1],
                                      params_mld_pike_night$SE[2],
                                      params_mld_pike_night$SE[3],
                                      params_mld_pike_night$SE[4],
                                      params_mld_pike_night$SE[5]),
                               ci_05 = c(params_mld_pike_night$CI_low[1],
                                         params_mld_pike_night$CI_low[2],
                                         params_mld_pike_night$CI_low[3],
                                         params_mld_pike_night$CI_low[4],
                                         params_mld_pike_night$CI_low[5]),
                               ci_95 = c(params_mld_pike_night$CI_high[1],
                                         params_mld_pike_night$CI_high[2],
                                         params_mld_pike_night$CI_high[3],
                                         params_mld_pike_night$CI_high[4],
                                         params_mld_pike_night$CI_high[5]),
                               t_value = c(params_mld_pike_night$t[1],
                                           params_mld_pike_night$t[2],
                                           params_mld_pike_night$t[3],
                                           params_mld_pike_night$t[4],
                                           params_mld_pike_night$t[5]),
                               p_value = c(params_mld_pike_night$p[1],
                                           params_mld_pike_night$p[2],
                                           params_mld_pike_night$p[3],
                                           params_mld_pike_night$p[4],
                                           params_mld_pike_night$p[5]))

#Random effects summary - Pike - Night ####
rdm_mdl_pike_night <- parameters::random_parameters(mdl_random_slo_int_night_pike)

mdl_pike_night_random <- tibble(species = "pike",
                                diel_period = "night",
                                description = c(rdm_mdl_pike_night$Description[1],
                                                rdm_mdl_pike_night$Description[2],
                                                rdm_mdl_pike_night$Description[3],
                                                rdm_mdl_pike_night$Description[4],
                                                rdm_mdl_pike_night$Description[5],
                                                rdm_mdl_pike_night$Description[6],
                                                rdm_mdl_pike_night$Description[7],
                                                rdm_mdl_pike_night$Description[8],
                                                rdm_mdl_pike_night$Description[9],
                                                rdm_mdl_pike_night$Description[10],
                                                rdm_mdl_pike_night$Description[11],
                                                rdm_mdl_pike_night$Description[12],
                                                rdm_mdl_pike_night$Description[13],
                                                rdm_mdl_pike_night$Description[14],
                                                rdm_mdl_pike_night$Description[15],
                                                rdm_mdl_pike_night$Description[16],
                                                rdm_mdl_pike_night$Description[17],
                                                rdm_mdl_pike_night$Description[18]),
                                component = c(rdm_mdl_pike_night$Component[1],
                                              rdm_mdl_pike_night$Component[2],
                                              rdm_mdl_pike_night$Component[3],
                                              rdm_mdl_pike_night$Component[4],
                                              rdm_mdl_pike_night$Component[5],
                                              rdm_mdl_pike_night$Component[6],
                                              rdm_mdl_pike_night$Component[7],
                                              rdm_mdl_pike_night$Component[8],
                                              rdm_mdl_pike_night$Component[9],
                                              rdm_mdl_pike_night$Component[10],
                                              rdm_mdl_pike_night$Component[11],
                                              rdm_mdl_pike_night$Component[12],
                                              rdm_mdl_pike_night$Component[13],
                                              rdm_mdl_pike_night$Component[14],
                                              rdm_mdl_pike_night$Component[15],
                                              rdm_mdl_pike_night$Component[16],
                                              rdm_mdl_pike_night$Component[17],
                                              rdm_mdl_pike_night$Component[18]),
                                type = c(rdm_mdl_pike_night$Type[1],
                                         rdm_mdl_pike_night$Type[2],
                                         rdm_mdl_pike_night$Type[3],
                                         rdm_mdl_pike_night$Type[4],
                                         rdm_mdl_pike_night$Type[5],
                                         rdm_mdl_pike_night$Type[6],
                                         rdm_mdl_pike_night$Type[7],
                                         rdm_mdl_pike_night$Type[8],
                                         rdm_mdl_pike_night$Type[9],
                                         rdm_mdl_pike_night$Type[10],
                                         rdm_mdl_pike_night$Type[11],
                                         rdm_mdl_pike_night$Type[12],
                                         rdm_mdl_pike_night$Type[13],
                                         rdm_mdl_pike_night$Type[14],
                                         rdm_mdl_pike_night$Type[15],
                                         rdm_mdl_pike_night$Type[16],
                                         rdm_mdl_pike_night$Type[17],
                                         rdm_mdl_pike_night$Type[18]),
                                term = c(rdm_mdl_pike_night$Term[1],
                                         rdm_mdl_pike_night$Term[2],
                                         rdm_mdl_pike_night$Term[3],
                                         rdm_mdl_pike_night$Term[4],
                                         rdm_mdl_pike_night$Term[5],
                                         rdm_mdl_pike_night$Term[6],
                                         rdm_mdl_pike_night$Term[7],
                                         rdm_mdl_pike_night$Term[8],
                                         rdm_mdl_pike_night$Term[9],
                                         rdm_mdl_pike_night$Term[10],
                                         rdm_mdl_pike_night$Term[11],
                                         rdm_mdl_pike_night$Term[12],
                                         rdm_mdl_pike_night$Term[13],
                                         rdm_mdl_pike_night$Term[14],
                                         rdm_mdl_pike_night$Term[15],
                                         rdm_mdl_pike_night$Term[16],
                                         rdm_mdl_pike_night$Term[17],
                                         rdm_mdl_pike_night$Term[18]),
                                value = c(rdm_mdl_pike_night$Value[1],
                                          rdm_mdl_pike_night$Value[2],
                                          rdm_mdl_pike_night$Value[3],
                                          rdm_mdl_pike_night$Value[4],
                                          rdm_mdl_pike_night$Value[5],
                                          rdm_mdl_pike_night$Value[6],
                                          rdm_mdl_pike_night$Value[7],
                                          rdm_mdl_pike_night$Value[8],
                                          rdm_mdl_pike_night$Value[9],
                                          rdm_mdl_pike_night$Value[10],
                                          rdm_mdl_pike_night$Value[11],
                                          rdm_mdl_pike_night$Value[12],
                                          rdm_mdl_pike_night$Value[13],
                                          rdm_mdl_pike_night$Value[14],
                                          rdm_mdl_pike_night$Value[15],
                                          rdm_mdl_pike_night$Value[16],
                                          rdm_mdl_pike_night$Value[17],
                                          rdm_mdl_pike_night$Value[18]))

#Percentage explained by each variable - Pike - Night ####
#Night
pike_mdl_night_percentage <- tibble(species = "pike",
                                    diel_period = "night",
                                    variable = c('Thermocline thickness (m)',
                                                 'Thermocline strength (ºC)',
                                                 'Thermocline depth (m)',
                                                 'Seiche strength'),
                                    percentage_explained = c(
                                      round(mdl_random_slo_int_pike_night_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_pike_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_pike_night_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_pike_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_pike_night_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_pike_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_pike_night_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_pike_night_aov$`Sum Sq`),2))
)

pike_night_mdl_summary <- list(fixed_effects = mdl_pike_night_fixed,
                               percentage_explained = pike_mdl_night_percentage,
                               random_effects = mdl_pike_night_random)

#Summary of the pike models - Pike####
pike_mld_summary <- list(fixed_effects = bind_rows(mdl_pike_day_fixed,
                                                   mdl_pike_night_fixed),
                         random_effects = bind_rows(mdl_pike_day_random,
                                                    mdl_pike_night_random),
                         percentage_explained = bind_rows(pike_mdl_day_percentage,
                                                          pike_mdl_night_percentage))

#Wels catfish ####
#Models run - Wels ####
#random intercept and random slope per fish id - Day (Seiche only)
tic('It took this amount of time to run the model')
mdl_random_slo_int_day_wels <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + 
                                 det_therm_strength + lake_therm_depth_smoothed_center + 
                                 det_therm_deviation_center + 
                                 (1 + lake_therm_thickness_smoothed + 
                                    det_therm_strength + 
                                    lake_therm_depth_smoothed_center +
                                    det_therm_deviation_center | fishid),
                               data = detections %>% 
                                 filter(species == "wels" & 
                                          diel_period == 'day' & 
                                          is_valid_seiche == TRUE), 
                               REML = T, 
                               lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()

#random intercept and random slope per fish id - Night (Seiche only)
tic('It took this amount of time to run the model')
mdl_random_slo_int_night_wels <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + 
                                   det_therm_strength + lake_therm_depth_smoothed_center + 
                                   det_therm_deviation_center + 
                                   (1 + lake_therm_thickness_smoothed + 
                                      det_therm_strength + 
                                      lake_therm_depth_smoothed_center + 
                                      det_therm_deviation_center | fishid),
                                 data = detections %>% 
                                   filter(species == "wels" & 
                                            diel_period == 'night' & 
                                            is_valid_seiche == TRUE), 
                                 REML = T, 
                                 lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()

#Summary of the model results - Wels ####
tab_model(mdl_random_slo_int_day_wels, 
          mdl_random_slo_int_night_wels, dv.labels = c("Day","Night"), 
          title='Wels catfish (Seiche only)',
          pred.labels = c('Intercept', 'Thermocline thickness (m)',
                          'Thermocline strength (ºC)',
                          'Thermocline depth (m)',
                          'Seiche strength'))#sigma² = withing group variance; tau00 = between group variance; ICC = proportion of variance explained for each group level

#Fixed effects summary - Wels - Day ####
params_mld_wels_day <- parameters::model_parameters(mdl_random_slo_int_day_wels)

mdl_wels_day_fixed <- tibble(species = "wels", 
                              diel_period = "day",
                              parameter = c(params_mld_wels_day$Parameter[1],
                                            params_mld_wels_day$Parameter[2],
                                            params_mld_wels_day$Parameter[3],
                                            params_mld_wels_day$Parameter[4],
                                            params_mld_wels_day$Parameter[5]),
                              coefficient = c(params_mld_wels_day$Coefficient[1],
                                              params_mld_wels_day$Coefficient[2],
                                              params_mld_wels_day$Coefficient[3],
                                              params_mld_wels_day$Coefficient[4],
                                              params_mld_wels_day$Coefficient[5]),
                              se = c(params_mld_wels_day$SE[1],
                                     params_mld_wels_day$SE[2],
                                     params_mld_wels_day$SE[3],
                                     params_mld_wels_day$SE[4],
                                     params_mld_wels_day$SE[5]),
                              ci_05 = c(params_mld_wels_day$CI_low[1],
                                        params_mld_wels_day$CI_low[2],
                                        params_mld_wels_day$CI_low[3],
                                        params_mld_wels_day$CI_low[4],
                                        params_mld_wels_day$CI_low[5]),
                              ci_95 = c(params_mld_wels_day$CI_high[1],
                                        params_mld_wels_day$CI_high[2],
                                        params_mld_wels_day$CI_high[3],
                                        params_mld_wels_day$CI_high[4],
                                        params_mld_wels_day$CI_high[5]),
                              t_value = c(params_mld_wels_day$t[1],
                                          params_mld_wels_day$t[2],
                                          params_mld_wels_day$t[3],
                                          params_mld_wels_day$t[4],
                                          params_mld_wels_day$t[5]),
                              p_value = c(params_mld_wels_day$p[1],
                                          params_mld_wels_day$p[2],
                                          params_mld_wels_day$p[3],
                                          params_mld_wels_day$p[4],
                                          params_mld_wels_day$p[5]))

#Random effects summary - Wels - Day ####
rdm_mdl_wels_day <- parameters::random_parameters(mdl_random_slo_int_day_wels)

mdl_wels_day_random <- tibble(species = "wels",
                               diel_period = "day",
                               description = c(rdm_mdl_wels_day$Description[1],
                                               rdm_mdl_wels_day$Description[2],
                                               rdm_mdl_wels_day$Description[3],
                                               rdm_mdl_wels_day$Description[4],
                                               rdm_mdl_wels_day$Description[5],
                                               rdm_mdl_wels_day$Description[6],
                                               rdm_mdl_wels_day$Description[7],
                                               rdm_mdl_wels_day$Description[8],
                                               rdm_mdl_wels_day$Description[9],
                                               rdm_mdl_wels_day$Description[10],
                                               rdm_mdl_wels_day$Description[11],
                                               rdm_mdl_wels_day$Description[12],
                                               rdm_mdl_wels_day$Description[13],
                                               rdm_mdl_wels_day$Description[14],
                                               rdm_mdl_wels_day$Description[15],
                                               rdm_mdl_wels_day$Description[16],
                                               rdm_mdl_wels_day$Description[17],
                                               rdm_mdl_wels_day$Description[18]),
                               component = c(rdm_mdl_wels_day$Component[1],
                                             rdm_mdl_wels_day$Component[2],
                                             rdm_mdl_wels_day$Component[3],
                                             rdm_mdl_wels_day$Component[4],
                                             rdm_mdl_wels_day$Component[5],
                                             rdm_mdl_wels_day$Component[6],
                                             rdm_mdl_wels_day$Component[7],
                                             rdm_mdl_wels_day$Component[8],
                                             rdm_mdl_wels_day$Component[9],
                                             rdm_mdl_wels_day$Component[10],
                                             rdm_mdl_wels_day$Component[11],
                                             rdm_mdl_wels_day$Component[12],
                                             rdm_mdl_wels_day$Component[13],
                                             rdm_mdl_wels_day$Component[14],
                                             rdm_mdl_wels_day$Component[15],
                                             rdm_mdl_wels_day$Component[16],
                                             rdm_mdl_wels_day$Component[17],
                                             rdm_mdl_wels_day$Component[18]),
                               type = c(rdm_mdl_wels_day$Type[1],
                                        rdm_mdl_wels_day$Type[2],
                                        rdm_mdl_wels_day$Type[3],
                                        rdm_mdl_wels_day$Type[4],
                                        rdm_mdl_wels_day$Type[5],
                                        rdm_mdl_wels_day$Type[6],
                                        rdm_mdl_wels_day$Type[7],
                                        rdm_mdl_wels_day$Type[8],
                                        rdm_mdl_wels_day$Type[9],
                                        rdm_mdl_wels_day$Type[10],
                                        rdm_mdl_wels_day$Type[11],
                                        rdm_mdl_wels_day$Type[12],
                                        rdm_mdl_wels_day$Type[13],
                                        rdm_mdl_wels_day$Type[14],
                                        rdm_mdl_wels_day$Type[15],
                                        rdm_mdl_wels_day$Type[16],
                                        rdm_mdl_wels_day$Type[17],
                                        rdm_mdl_wels_day$Type[18]),
                               term = c(rdm_mdl_wels_day$Term[1],
                                        rdm_mdl_wels_day$Term[2],
                                        rdm_mdl_wels_day$Term[3],
                                        rdm_mdl_wels_day$Term[4],
                                        rdm_mdl_wels_day$Term[5],
                                        rdm_mdl_wels_day$Term[6],
                                        rdm_mdl_wels_day$Term[7],
                                        rdm_mdl_wels_day$Term[8],
                                        rdm_mdl_wels_day$Term[9],
                                        rdm_mdl_wels_day$Term[10],
                                        rdm_mdl_wels_day$Term[11],
                                        rdm_mdl_wels_day$Term[12],
                                        rdm_mdl_wels_day$Term[13],
                                        rdm_mdl_wels_day$Term[14],
                                        rdm_mdl_wels_day$Term[15],
                                        rdm_mdl_wels_day$Term[16],
                                        rdm_mdl_wels_day$Term[17],
                                        rdm_mdl_wels_day$Term[18]),
                               value = c(rdm_mdl_wels_day$Value[1],
                                         rdm_mdl_wels_day$Value[2],
                                         rdm_mdl_wels_day$Value[3],
                                         rdm_mdl_wels_day$Value[4],
                                         rdm_mdl_wels_day$Value[5],
                                         rdm_mdl_wels_day$Value[6],
                                         rdm_mdl_wels_day$Value[7],
                                         rdm_mdl_wels_day$Value[8],
                                         rdm_mdl_wels_day$Value[9],
                                         rdm_mdl_wels_day$Value[10],
                                         rdm_mdl_wels_day$Value[11],
                                         rdm_mdl_wels_day$Value[12],
                                         rdm_mdl_wels_day$Value[13],
                                         rdm_mdl_wels_day$Value[14],
                                         rdm_mdl_wels_day$Value[15],
                                         rdm_mdl_wels_day$Value[16],
                                         rdm_mdl_wels_day$Value[17],
                                         rdm_mdl_wels_day$Value[18]))

#Percentage explained by each variable - Wels - Day ####
#Sum of squares to give the percentage explained by each variable in the model
#Day
mdl_random_slo_int_wels_day_aov <- anova(mdl_random_slo_int_day_wels)
wels_mdl_day_percentage <- tibble(species = "wels",
                                   diel_period = "day",
                                   variable = c('Thermocline thickness (m)',
                                                'Thermocline strength (ºC)',
                                                'Thermocline depth (m)',
                                                'Seiche strength'),
                                   percentage_explained = c(
                                     round(mdl_random_slo_int_wels_day_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_wels_day_aov$`Sum Sq`),2),
                                     round(mdl_random_slo_int_wels_day_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_wels_day_aov$`Sum Sq`),2),
                                     round(mdl_random_slo_int_wels_day_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_wels_day_aov$`Sum Sq`),2),
                                     round(mdl_random_slo_int_wels_day_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_wels_day_aov$`Sum Sq`),2)
                                   )
)

wels_day_mdl_summary <- list(fixed_effects = mdl_wels_day_fixed,
                              percentage_explained = wels_mdl_day_percentage,
                              random_effects = mdl_wels_day_random)

#Night
mdl_random_slo_int_wels_night_aov <-  anova(mdl_random_slo_int_night_wels)
wels_mdl_night_percentage <- tibble(species = "wels",
                                     diel_period = "night",
                                     variable = c('Thermocline thickness (m)',
                                                  'Thermocline strength (ºC)',
                                                  'Thermocline depth (m)',
                                                  'Seiche strength'),
                                     percentage_explained = c(
                                       round(mdl_random_slo_int_wels_night_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_wels_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_wels_night_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_wels_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_wels_night_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_wels_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_wels_night_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_wels_night_aov$`Sum Sq`),2))
)

#wels (day & night)
wels_mld_percentages <- bind_rows(wels_mdl_day_percentage, wels_mdl_night_percentage)

#Fixed effects summary - Wels - Night ####
params_mld_wels_night <- parameters::model_parameters(mdl_random_slo_int_night_wels)

mdl_wels_night_fixed <- tibble(species = "wels", 
                                diel_period = "night",
                                parameter = c(params_mld_wels_night$Parameter[1],
                                              params_mld_wels_night$Parameter[2],
                                              params_mld_wels_night$Parameter[3],
                                              params_mld_wels_night$Parameter[4],
                                              params_mld_wels_night$Parameter[5]),
                                coefficient = c(params_mld_wels_night$Coefficient[1],
                                                params_mld_wels_night$Coefficient[2],
                                                params_mld_wels_night$Coefficient[3],
                                                params_mld_wels_night$Coefficient[4],
                                                params_mld_wels_night$Coefficient[5]),
                                se = c(params_mld_wels_night$SE[1],
                                       params_mld_wels_night$SE[2],
                                       params_mld_wels_night$SE[3],
                                       params_mld_wels_night$SE[4],
                                       params_mld_wels_night$SE[5]),
                                ci_05 = c(params_mld_wels_night$CI_low[1],
                                          params_mld_wels_night$CI_low[2],
                                          params_mld_wels_night$CI_low[3],
                                          params_mld_wels_night$CI_low[4],
                                          params_mld_wels_night$CI_low[5]),
                                ci_95 = c(params_mld_wels_night$CI_high[1],
                                          params_mld_wels_night$CI_high[2],
                                          params_mld_wels_night$CI_high[3],
                                          params_mld_wels_night$CI_high[4],
                                          params_mld_wels_night$CI_high[5]),
                                t_value = c(params_mld_wels_night$t[1],
                                            params_mld_wels_night$t[2],
                                            params_mld_wels_night$t[3],
                                            params_mld_wels_night$t[4],
                                            params_mld_wels_night$t[5]),
                                p_value = c(params_mld_wels_night$p[1],
                                            params_mld_wels_night$p[2],
                                            params_mld_wels_night$p[3],
                                            params_mld_wels_night$p[4],
                                            params_mld_wels_night$p[5]))

#Random effects summary - Wels - Night ####
rdm_mdl_wels_night <- parameters::random_parameters(mdl_random_slo_int_night_wels)

mdl_wels_night_random <- tibble(species = "wels",
                                 diel_period = "night",
                                 description = c(rdm_mdl_wels_night$Description[1],
                                                 rdm_mdl_wels_night$Description[2],
                                                 rdm_mdl_wels_night$Description[3],
                                                 rdm_mdl_wels_night$Description[4],
                                                 rdm_mdl_wels_night$Description[5],
                                                 rdm_mdl_wels_night$Description[6],
                                                 rdm_mdl_wels_night$Description[7],
                                                 rdm_mdl_wels_night$Description[8],
                                                 rdm_mdl_wels_night$Description[9],
                                                 rdm_mdl_wels_night$Description[10],
                                                 rdm_mdl_wels_night$Description[11],
                                                 rdm_mdl_wels_night$Description[12],
                                                 rdm_mdl_wels_night$Description[13],
                                                 rdm_mdl_wels_night$Description[14],
                                                 rdm_mdl_wels_night$Description[15],
                                                 rdm_mdl_wels_night$Description[16],
                                                 rdm_mdl_wels_night$Description[17],
                                                 rdm_mdl_wels_night$Description[18]),
                                 component = c(rdm_mdl_wels_night$Component[1],
                                               rdm_mdl_wels_night$Component[2],
                                               rdm_mdl_wels_night$Component[3],
                                               rdm_mdl_wels_night$Component[4],
                                               rdm_mdl_wels_night$Component[5],
                                               rdm_mdl_wels_night$Component[6],
                                               rdm_mdl_wels_night$Component[7],
                                               rdm_mdl_wels_night$Component[8],
                                               rdm_mdl_wels_night$Component[9],
                                               rdm_mdl_wels_night$Component[10],
                                               rdm_mdl_wels_night$Component[11],
                                               rdm_mdl_wels_night$Component[12],
                                               rdm_mdl_wels_night$Component[13],
                                               rdm_mdl_wels_night$Component[14],
                                               rdm_mdl_wels_night$Component[15],
                                               rdm_mdl_wels_night$Component[16],
                                               rdm_mdl_wels_night$Component[17],
                                               rdm_mdl_wels_night$Component[18]),
                                 type = c(rdm_mdl_wels_night$Type[1],
                                          rdm_mdl_wels_night$Type[2],
                                          rdm_mdl_wels_night$Type[3],
                                          rdm_mdl_wels_night$Type[4],
                                          rdm_mdl_wels_night$Type[5],
                                          rdm_mdl_wels_night$Type[6],
                                          rdm_mdl_wels_night$Type[7],
                                          rdm_mdl_wels_night$Type[8],
                                          rdm_mdl_wels_night$Type[9],
                                          rdm_mdl_wels_night$Type[10],
                                          rdm_mdl_wels_night$Type[11],
                                          rdm_mdl_wels_night$Type[12],
                                          rdm_mdl_wels_night$Type[13],
                                          rdm_mdl_wels_night$Type[14],
                                          rdm_mdl_wels_night$Type[15],
                                          rdm_mdl_wels_night$Type[16],
                                          rdm_mdl_wels_night$Type[17],
                                          rdm_mdl_wels_night$Type[18]),
                                 term = c(rdm_mdl_wels_night$Term[1],
                                          rdm_mdl_wels_night$Term[2],
                                          rdm_mdl_wels_night$Term[3],
                                          rdm_mdl_wels_night$Term[4],
                                          rdm_mdl_wels_night$Term[5],
                                          rdm_mdl_wels_night$Term[6],
                                          rdm_mdl_wels_night$Term[7],
                                          rdm_mdl_wels_night$Term[8],
                                          rdm_mdl_wels_night$Term[9],
                                          rdm_mdl_wels_night$Term[10],
                                          rdm_mdl_wels_night$Term[11],
                                          rdm_mdl_wels_night$Term[12],
                                          rdm_mdl_wels_night$Term[13],
                                          rdm_mdl_wels_night$Term[14],
                                          rdm_mdl_wels_night$Term[15],
                                          rdm_mdl_wels_night$Term[16],
                                          rdm_mdl_wels_night$Term[17],
                                          rdm_mdl_wels_night$Term[18]),
                                 value = c(rdm_mdl_wels_night$Value[1],
                                           rdm_mdl_wels_night$Value[2],
                                           rdm_mdl_wels_night$Value[3],
                                           rdm_mdl_wels_night$Value[4],
                                           rdm_mdl_wels_night$Value[5],
                                           rdm_mdl_wels_night$Value[6],
                                           rdm_mdl_wels_night$Value[7],
                                           rdm_mdl_wels_night$Value[8],
                                           rdm_mdl_wels_night$Value[9],
                                           rdm_mdl_wels_night$Value[10],
                                           rdm_mdl_wels_night$Value[11],
                                           rdm_mdl_wels_night$Value[12],
                                           rdm_mdl_wels_night$Value[13],
                                           rdm_mdl_wels_night$Value[14],
                                           rdm_mdl_wels_night$Value[15],
                                           rdm_mdl_wels_night$Value[16],
                                           rdm_mdl_wels_night$Value[17],
                                           rdm_mdl_wels_night$Value[18]))

#Percentage explained by each variable - Wels - Night ####
#Night
wels_mdl_night_percentage <- tibble(species = "wels",
                                     diel_period = "night",
                                     variable = c('Thermocline thickness (m)',
                                                  'Thermocline strength (ºC)',
                                                  'Thermocline depth (m)',
                                                  'Seiche strength'),
                                     percentage_explained = c(
                                       round(mdl_random_slo_int_wels_night_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_wels_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_wels_night_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_wels_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_wels_night_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_wels_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_wels_night_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_wels_night_aov$`Sum Sq`),2))
)

wels_night_mdl_summary <- list(fixed_effects = mdl_wels_night_fixed,
                                percentage_explained = wels_mdl_night_percentage,
                                random_effects = mdl_wels_night_random)

#Summary of the wels models - Wels####
wels_mld_summary <- list(fixed_effects = bind_rows(mdl_wels_day_fixed,
                                                    mdl_wels_night_fixed),
                          random_effects = bind_rows(mdl_wels_day_random,
                                                     mdl_wels_night_random),
                          percentage_explained = bind_rows(wels_mdl_day_percentage,
                                                           wels_mdl_night_percentage))
#Tench ####
#Models run - Tench #####
#random intercept and random slope per fish id - Day (Seiche only)
tic('It took this amount of time to run the model')
mdl_random_slo_int_day_tench <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + 
                                      det_therm_strength + lake_therm_depth_smoothed_center + 
                                      det_therm_deviation_center + 
                                      (1 + lake_therm_thickness_smoothed + 
                                         det_therm_strength + 
                                         lake_therm_depth_smoothed_center +
                                         det_therm_deviation_center | fishid),
                                    data = detections %>% 
                                      filter(species == "tench" & 
                                               diel_period == 'day' & 
                                               is_valid_seiche == TRUE), 
                                    REML = T, 
                                    lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()

#random intercept and random slope per fish id - Night (Seiche only)
tic('It took this amount of time to run the model')
mdl_random_slo_int_night_tench <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + 
                                        det_therm_strength + lake_therm_depth_smoothed_center + 
                                        det_therm_deviation_center + 
                                        (1 + lake_therm_thickness_smoothed + 
                                           det_therm_strength + 
                                           lake_therm_depth_smoothed_center + 
                                           det_therm_deviation_center | fishid),
                                      data = detections %>% 
                                        filter(species == "tench" & 
                                                 diel_period == 'night' & 
                                                 is_valid_seiche == TRUE), 
                                      REML = T, 
                                      lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()

#Summary of the model results - Tench #####
tab_model(mdl_random_slo_int_day_tench, 
          mdl_random_slo_int_night_tench, dv.labels = c("Day","Night"), 
          title='tench (Seiche only)',
          pred.labels = c('Intercept', 'Thermocline thickness (m)',
                          'Thermocline strength (ºC)',
                          'Thermocline depth (m)',
                          'Seiche strength'))#sigma² = withing group variance; tau00 = between group variance; ICC = proportion of variance explained for each group level

#Fixed effects summary - Tench - Day ####
params_mld_tench_day <- parameters::model_parameters(mdl_random_slo_int_day_tench)

mdl_tench_day_fixed <- tibble(species = "tench", 
                             diel_period = "day",
                             parameter = c(params_mld_tench_day$Parameter[1],
                                           params_mld_tench_day$Parameter[2],
                                           params_mld_tench_day$Parameter[3],
                                           params_mld_tench_day$Parameter[4],
                                           params_mld_tench_day$Parameter[5]),
                             coefficient = c(params_mld_tench_day$Coefficient[1],
                                             params_mld_tench_day$Coefficient[2],
                                             params_mld_tench_day$Coefficient[3],
                                             params_mld_tench_day$Coefficient[4],
                                             params_mld_tench_day$Coefficient[5]),
                             se = c(params_mld_tench_day$SE[1],
                                    params_mld_tench_day$SE[2],
                                    params_mld_tench_day$SE[3],
                                    params_mld_tench_day$SE[4],
                                    params_mld_tench_day$SE[5]),
                             ci_05 = c(params_mld_tench_day$CI_low[1],
                                       params_mld_tench_day$CI_low[2],
                                       params_mld_tench_day$CI_low[3],
                                       params_mld_tench_day$CI_low[4],
                                       params_mld_tench_day$CI_low[5]),
                             ci_95 = c(params_mld_tench_day$CI_high[1],
                                       params_mld_tench_day$CI_high[2],
                                       params_mld_tench_day$CI_high[3],
                                       params_mld_tench_day$CI_high[4],
                                       params_mld_tench_day$CI_high[5]),
                             t_value = c(params_mld_tench_day$t[1],
                                         params_mld_tench_day$t[2],
                                         params_mld_tench_day$t[3],
                                         params_mld_tench_day$t[4],
                                         params_mld_tench_day$t[5]),
                             p_value = c(params_mld_tench_day$p[1],
                                         params_mld_tench_day$p[2],
                                         params_mld_tench_day$p[3],
                                         params_mld_tench_day$p[4],
                                         params_mld_tench_day$p[5]))

#Random effects summary - Tench - Day ####
rdm_mdl_tench_day <- parameters::random_parameters(mdl_random_slo_int_day_tench)

mdl_tench_day_random <- tibble(species = "tench",
                              diel_period = "day",
                              description = c(rdm_mdl_tench_day$Description[1],
                                              rdm_mdl_tench_day$Description[2],
                                              rdm_mdl_tench_day$Description[3],
                                              rdm_mdl_tench_day$Description[4],
                                              rdm_mdl_tench_day$Description[5],
                                              rdm_mdl_tench_day$Description[6],
                                              rdm_mdl_tench_day$Description[7],
                                              rdm_mdl_tench_day$Description[8],
                                              rdm_mdl_tench_day$Description[9],
                                              rdm_mdl_tench_day$Description[10],
                                              rdm_mdl_tench_day$Description[11],
                                              rdm_mdl_tench_day$Description[12],
                                              rdm_mdl_tench_day$Description[13],
                                              rdm_mdl_tench_day$Description[14],
                                              rdm_mdl_tench_day$Description[15],
                                              rdm_mdl_tench_day$Description[16],
                                              rdm_mdl_tench_day$Description[17],
                                              rdm_mdl_tench_day$Description[18]),
                              component = c(rdm_mdl_tench_day$Component[1],
                                            rdm_mdl_tench_day$Component[2],
                                            rdm_mdl_tench_day$Component[3],
                                            rdm_mdl_tench_day$Component[4],
                                            rdm_mdl_tench_day$Component[5],
                                            rdm_mdl_tench_day$Component[6],
                                            rdm_mdl_tench_day$Component[7],
                                            rdm_mdl_tench_day$Component[8],
                                            rdm_mdl_tench_day$Component[9],
                                            rdm_mdl_tench_day$Component[10],
                                            rdm_mdl_tench_day$Component[11],
                                            rdm_mdl_tench_day$Component[12],
                                            rdm_mdl_tench_day$Component[13],
                                            rdm_mdl_tench_day$Component[14],
                                            rdm_mdl_tench_day$Component[15],
                                            rdm_mdl_tench_day$Component[16],
                                            rdm_mdl_tench_day$Component[17],
                                            rdm_mdl_tench_day$Component[18]),
                              type = c(rdm_mdl_tench_day$Type[1],
                                       rdm_mdl_tench_day$Type[2],
                                       rdm_mdl_tench_day$Type[3],
                                       rdm_mdl_tench_day$Type[4],
                                       rdm_mdl_tench_day$Type[5],
                                       rdm_mdl_tench_day$Type[6],
                                       rdm_mdl_tench_day$Type[7],
                                       rdm_mdl_tench_day$Type[8],
                                       rdm_mdl_tench_day$Type[9],
                                       rdm_mdl_tench_day$Type[10],
                                       rdm_mdl_tench_day$Type[11],
                                       rdm_mdl_tench_day$Type[12],
                                       rdm_mdl_tench_day$Type[13],
                                       rdm_mdl_tench_day$Type[14],
                                       rdm_mdl_tench_day$Type[15],
                                       rdm_mdl_tench_day$Type[16],
                                       rdm_mdl_tench_day$Type[17],
                                       rdm_mdl_tench_day$Type[18]),
                              term = c(rdm_mdl_tench_day$Term[1],
                                       rdm_mdl_tench_day$Term[2],
                                       rdm_mdl_tench_day$Term[3],
                                       rdm_mdl_tench_day$Term[4],
                                       rdm_mdl_tench_day$Term[5],
                                       rdm_mdl_tench_day$Term[6],
                                       rdm_mdl_tench_day$Term[7],
                                       rdm_mdl_tench_day$Term[8],
                                       rdm_mdl_tench_day$Term[9],
                                       rdm_mdl_tench_day$Term[10],
                                       rdm_mdl_tench_day$Term[11],
                                       rdm_mdl_tench_day$Term[12],
                                       rdm_mdl_tench_day$Term[13],
                                       rdm_mdl_tench_day$Term[14],
                                       rdm_mdl_tench_day$Term[15],
                                       rdm_mdl_tench_day$Term[16],
                                       rdm_mdl_tench_day$Term[17],
                                       rdm_mdl_tench_day$Term[18]),
                              value = c(rdm_mdl_tench_day$Value[1],
                                        rdm_mdl_tench_day$Value[2],
                                        rdm_mdl_tench_day$Value[3],
                                        rdm_mdl_tench_day$Value[4],
                                        rdm_mdl_tench_day$Value[5],
                                        rdm_mdl_tench_day$Value[6],
                                        rdm_mdl_tench_day$Value[7],
                                        rdm_mdl_tench_day$Value[8],
                                        rdm_mdl_tench_day$Value[9],
                                        rdm_mdl_tench_day$Value[10],
                                        rdm_mdl_tench_day$Value[11],
                                        rdm_mdl_tench_day$Value[12],
                                        rdm_mdl_tench_day$Value[13],
                                        rdm_mdl_tench_day$Value[14],
                                        rdm_mdl_tench_day$Value[15],
                                        rdm_mdl_tench_day$Value[16],
                                        rdm_mdl_tench_day$Value[17],
                                        rdm_mdl_tench_day$Value[18]))

#Percentage explained by each variable - Tench - Day ####
#Sum of squares to give the percentage explained by each variable in the model
#Day
mdl_random_slo_int_tench_day_aov <- anova(mdl_random_slo_int_day_tench)
tench_mdl_day_percentage <- tibble(species = "tench",
                                  diel_period = "day",
                                  variable = c('Thermocline thickness (m)',
                                               'Thermocline strength (ºC)',
                                               'Thermocline depth (m)',
                                               'Seiche strength'),
                                  percentage_explained = c(
                                    round(mdl_random_slo_int_tench_day_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_tench_day_aov$`Sum Sq`),2),
                                    round(mdl_random_slo_int_tench_day_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_tench_day_aov$`Sum Sq`),2),
                                    round(mdl_random_slo_int_tench_day_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_tench_day_aov$`Sum Sq`),2),
                                    round(mdl_random_slo_int_tench_day_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_tench_day_aov$`Sum Sq`),2)
                                  )
)

tench_day_mdl_summary <- list(fixed_effects = mdl_tench_day_fixed,
                             percentage_explained = tench_mdl_day_percentage,
                             random_effects = mdl_tench_day_random)

#Night
mdl_random_slo_int_tench_night_aov <-  anova(mdl_random_slo_int_night_tench)
tench_mdl_night_percentage <- tibble(species = "tench",
                                    diel_period = "night",
                                    variable = c('Thermocline thickness (m)',
                                                 'Thermocline strength (ºC)',
                                                 'Thermocline depth (m)',
                                                 'Seiche strength'),
                                    percentage_explained = c(
                                      round(mdl_random_slo_int_tench_night_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_tench_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_tench_night_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_tench_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_tench_night_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_tench_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_tench_night_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_tench_night_aov$`Sum Sq`),2))
)

#tench (day & night)
tench_mld_percentages <- bind_rows(tench_mdl_day_percentage, tench_mdl_night_percentage)

#Fixed effects summary - Tench - Night ####
params_mld_tench_night <- parameters::model_parameters(mdl_random_slo_int_night_tench)

mdl_tench_night_fixed <- tibble(species = "tench", 
                               diel_period = "night",
                               parameter = c(params_mld_tench_night$Parameter[1],
                                             params_mld_tench_night$Parameter[2],
                                             params_mld_tench_night$Parameter[3],
                                             params_mld_tench_night$Parameter[4],
                                             params_mld_tench_night$Parameter[5]),
                               coefficient = c(params_mld_tench_night$Coefficient[1],
                                               params_mld_tench_night$Coefficient[2],
                                               params_mld_tench_night$Coefficient[3],
                                               params_mld_tench_night$Coefficient[4],
                                               params_mld_tench_night$Coefficient[5]),
                               se = c(params_mld_tench_night$SE[1],
                                      params_mld_tench_night$SE[2],
                                      params_mld_tench_night$SE[3],
                                      params_mld_tench_night$SE[4],
                                      params_mld_tench_night$SE[5]),
                               ci_05 = c(params_mld_tench_night$CI_low[1],
                                         params_mld_tench_night$CI_low[2],
                                         params_mld_tench_night$CI_low[3],
                                         params_mld_tench_night$CI_low[4],
                                         params_mld_tench_night$CI_low[5]),
                               ci_95 = c(params_mld_tench_night$CI_high[1],
                                         params_mld_tench_night$CI_high[2],
                                         params_mld_tench_night$CI_high[3],
                                         params_mld_tench_night$CI_high[4],
                                         params_mld_tench_night$CI_high[5]),
                               t_value = c(params_mld_tench_night$t[1],
                                           params_mld_tench_night$t[2],
                                           params_mld_tench_night$t[3],
                                           params_mld_tench_night$t[4],
                                           params_mld_tench_night$t[5]),
                               p_value = c(params_mld_tench_night$p[1],
                                           params_mld_tench_night$p[2],
                                           params_mld_tench_night$p[3],
                                           params_mld_tench_night$p[4],
                                           params_mld_tench_night$p[5]))

#Random effects summary - Tench - Night ####
rdm_mdl_tench_night <- parameters::random_parameters(mdl_random_slo_int_night_tench)

mdl_tench_night_random <- tibble(species = "tench",
                                diel_period = "night",
                                description = c(rdm_mdl_tench_night$Description[1],
                                                rdm_mdl_tench_night$Description[2],
                                                rdm_mdl_tench_night$Description[3],
                                                rdm_mdl_tench_night$Description[4],
                                                rdm_mdl_tench_night$Description[5],
                                                rdm_mdl_tench_night$Description[6],
                                                rdm_mdl_tench_night$Description[7],
                                                rdm_mdl_tench_night$Description[8],
                                                rdm_mdl_tench_night$Description[9],
                                                rdm_mdl_tench_night$Description[10],
                                                rdm_mdl_tench_night$Description[11],
                                                rdm_mdl_tench_night$Description[12],
                                                rdm_mdl_tench_night$Description[13],
                                                rdm_mdl_tench_night$Description[14],
                                                rdm_mdl_tench_night$Description[15],
                                                rdm_mdl_tench_night$Description[16],
                                                rdm_mdl_tench_night$Description[17],
                                                rdm_mdl_tench_night$Description[18]),
                                component = c(rdm_mdl_tench_night$Component[1],
                                              rdm_mdl_tench_night$Component[2],
                                              rdm_mdl_tench_night$Component[3],
                                              rdm_mdl_tench_night$Component[4],
                                              rdm_mdl_tench_night$Component[5],
                                              rdm_mdl_tench_night$Component[6],
                                              rdm_mdl_tench_night$Component[7],
                                              rdm_mdl_tench_night$Component[8],
                                              rdm_mdl_tench_night$Component[9],
                                              rdm_mdl_tench_night$Component[10],
                                              rdm_mdl_tench_night$Component[11],
                                              rdm_mdl_tench_night$Component[12],
                                              rdm_mdl_tench_night$Component[13],
                                              rdm_mdl_tench_night$Component[14],
                                              rdm_mdl_tench_night$Component[15],
                                              rdm_mdl_tench_night$Component[16],
                                              rdm_mdl_tench_night$Component[17],
                                              rdm_mdl_tench_night$Component[18]),
                                type = c(rdm_mdl_tench_night$Type[1],
                                         rdm_mdl_tench_night$Type[2],
                                         rdm_mdl_tench_night$Type[3],
                                         rdm_mdl_tench_night$Type[4],
                                         rdm_mdl_tench_night$Type[5],
                                         rdm_mdl_tench_night$Type[6],
                                         rdm_mdl_tench_night$Type[7],
                                         rdm_mdl_tench_night$Type[8],
                                         rdm_mdl_tench_night$Type[9],
                                         rdm_mdl_tench_night$Type[10],
                                         rdm_mdl_tench_night$Type[11],
                                         rdm_mdl_tench_night$Type[12],
                                         rdm_mdl_tench_night$Type[13],
                                         rdm_mdl_tench_night$Type[14],
                                         rdm_mdl_tench_night$Type[15],
                                         rdm_mdl_tench_night$Type[16],
                                         rdm_mdl_tench_night$Type[17],
                                         rdm_mdl_tench_night$Type[18]),
                                term = c(rdm_mdl_tench_night$Term[1],
                                         rdm_mdl_tench_night$Term[2],
                                         rdm_mdl_tench_night$Term[3],
                                         rdm_mdl_tench_night$Term[4],
                                         rdm_mdl_tench_night$Term[5],
                                         rdm_mdl_tench_night$Term[6],
                                         rdm_mdl_tench_night$Term[7],
                                         rdm_mdl_tench_night$Term[8],
                                         rdm_mdl_tench_night$Term[9],
                                         rdm_mdl_tench_night$Term[10],
                                         rdm_mdl_tench_night$Term[11],
                                         rdm_mdl_tench_night$Term[12],
                                         rdm_mdl_tench_night$Term[13],
                                         rdm_mdl_tench_night$Term[14],
                                         rdm_mdl_tench_night$Term[15],
                                         rdm_mdl_tench_night$Term[16],
                                         rdm_mdl_tench_night$Term[17],
                                         rdm_mdl_tench_night$Term[18]),
                                value = c(rdm_mdl_tench_night$Value[1],
                                          rdm_mdl_tench_night$Value[2],
                                          rdm_mdl_tench_night$Value[3],
                                          rdm_mdl_tench_night$Value[4],
                                          rdm_mdl_tench_night$Value[5],
                                          rdm_mdl_tench_night$Value[6],
                                          rdm_mdl_tench_night$Value[7],
                                          rdm_mdl_tench_night$Value[8],
                                          rdm_mdl_tench_night$Value[9],
                                          rdm_mdl_tench_night$Value[10],
                                          rdm_mdl_tench_night$Value[11],
                                          rdm_mdl_tench_night$Value[12],
                                          rdm_mdl_tench_night$Value[13],
                                          rdm_mdl_tench_night$Value[14],
                                          rdm_mdl_tench_night$Value[15],
                                          rdm_mdl_tench_night$Value[16],
                                          rdm_mdl_tench_night$Value[17],
                                          rdm_mdl_tench_night$Value[18]))
#Percentage explained by each variable - Tench - Night ####
#Night
tench_mdl_night_percentage <- tibble(species = "tench",
                                    diel_period = "night",
                                    variable = c('Thermocline thickness (m)',
                                                 'Thermocline strength (ºC)',
                                                 'Thermocline depth (m)',
                                                 'Seiche strength'),
                                    percentage_explained = c(
                                      round(mdl_random_slo_int_tench_night_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_tench_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_tench_night_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_tench_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_tench_night_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_tench_night_aov$`Sum Sq`),2),
                                      round(mdl_random_slo_int_tench_night_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_tench_night_aov$`Sum Sq`),2))
)


tench_night_mdl_summary <- list(fixed_effects = mdl_tench_night_fixed,
                               percentage_explained = tench_mdl_night_percentage,
                               random_effects = mdl_tench_night_random)

#Summary of the tench models - Tench####
tench_mld_summary <- list(fixed_effects = bind_rows(mdl_tench_day_fixed,
                                                   mdl_tench_night_fixed),
                         random_effects = bind_rows(mdl_tench_day_random,
                                                    mdl_tench_night_random),
                         percentage_explained = bind_rows(tench_mdl_day_percentage,
                                                          tench_mdl_night_percentage))

#Rudd ####
#Models run - Rudd ####
#random intercept and random slope per fish id - Day (Seiche only)
tic('It took this amount of time to run the model')
mdl_random_slo_int_day_rudd <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + 
                                       det_therm_strength + lake_therm_depth_smoothed_center + 
                                       det_therm_deviation_center + 
                                       (1 + lake_therm_thickness_smoothed + 
                                          det_therm_strength + 
                                          lake_therm_depth_smoothed_center +
                                          det_therm_deviation_center | fishid),
                                     data = detections %>% 
                                       filter(species == "rudd" & 
                                                diel_period == 'day' & 
                                                is_valid_seiche == TRUE), 
                                     REML = T, 
                                     lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()

#random intercept and random slope per fish id - Night (Seiche only)
tic('It took this amount of time to run the model')
mdl_random_slo_int_night_rudd <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + 
                                         det_therm_strength + lake_therm_depth_smoothed_center + 
                                         det_therm_deviation_center + 
                                         (1 + lake_therm_thickness_smoothed + 
                                            det_therm_strength + 
                                            lake_therm_depth_smoothed_center + 
                                            det_therm_deviation_center | fishid),
                                       data = detections %>% 
                                         filter(species == "rudd" & 
                                                  diel_period == 'night' & 
                                                  is_valid_seiche == TRUE), 
                                       REML = T, 
                                       lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()

#Summary of the model results - Rudd #####
tab_model(mdl_random_slo_int_day_rudd, 
          mdl_random_slo_int_night_rudd, dv.labels = c("Day","Night"), 
          title='Rudd (Seiche only)',
          pred.labels = c('Intercept', 'Thermocline thickness (m)',
                          'Thermocline strength (ºC)',
                          'Thermocline depth (m)',
                          'Seiche strength'))#sigma² = withing group variance; tau00 = between group variance; ICC = proportion of variance explained for each group level

#Fixed effects summary - Rudd - Day ####
params_mld_rudd_day <- parameters::model_parameters(mdl_random_slo_int_day_rudd)

mdl_rudd_day_fixed <- tibble(species = "rudd", 
                              diel_period = "day",
                              parameter = c(params_mld_rudd_day$Parameter[1],
                                            params_mld_rudd_day$Parameter[2],
                                            params_mld_rudd_day$Parameter[3],
                                            params_mld_rudd_day$Parameter[4],
                                            params_mld_rudd_day$Parameter[5]),
                              coefficient = c(params_mld_rudd_day$Coefficient[1],
                                              params_mld_rudd_day$Coefficient[2],
                                              params_mld_rudd_day$Coefficient[3],
                                              params_mld_rudd_day$Coefficient[4],
                                              params_mld_rudd_day$Coefficient[5]),
                              se = c(params_mld_rudd_day$SE[1],
                                     params_mld_rudd_day$SE[2],
                                     params_mld_rudd_day$SE[3],
                                     params_mld_rudd_day$SE[4],
                                     params_mld_rudd_day$SE[5]),
                              ci_05 = c(params_mld_rudd_day$CI_low[1],
                                        params_mld_rudd_day$CI_low[2],
                                        params_mld_rudd_day$CI_low[3],
                                        params_mld_rudd_day$CI_low[4],
                                        params_mld_rudd_day$CI_low[5]),
                              ci_95 = c(params_mld_rudd_day$CI_high[1],
                                        params_mld_rudd_day$CI_high[2],
                                        params_mld_rudd_day$CI_high[3],
                                        params_mld_rudd_day$CI_high[4],
                                        params_mld_rudd_day$CI_high[5]),
                              t_value = c(params_mld_rudd_day$t[1],
                                          params_mld_rudd_day$t[2],
                                          params_mld_rudd_day$t[3],
                                          params_mld_rudd_day$t[4],
                                          params_mld_rudd_day$t[5]),
                              p_value = c(params_mld_rudd_day$p[1],
                                          params_mld_rudd_day$p[2],
                                          params_mld_rudd_day$p[3],
                                          params_mld_rudd_day$p[4],
                                          params_mld_rudd_day$p[5]))

#Random effects summary - Rudd - Day ####
rdm_mdl_rudd_day <- parameters::random_parameters(mdl_random_slo_int_day_rudd)

mdl_rudd_day_random <- tibble(species = "rudd",
                               diel_period = "day",
                               description = c(rdm_mdl_rudd_day$Description[1],
                                               rdm_mdl_rudd_day$Description[2],
                                               rdm_mdl_rudd_day$Description[3],
                                               rdm_mdl_rudd_day$Description[4],
                                               rdm_mdl_rudd_day$Description[5],
                                               rdm_mdl_rudd_day$Description[6],
                                               rdm_mdl_rudd_day$Description[7],
                                               rdm_mdl_rudd_day$Description[8],
                                               rdm_mdl_rudd_day$Description[9],
                                               rdm_mdl_rudd_day$Description[10],
                                               rdm_mdl_rudd_day$Description[11],
                                               rdm_mdl_rudd_day$Description[12],
                                               rdm_mdl_rudd_day$Description[13],
                                               rdm_mdl_rudd_day$Description[14],
                                               rdm_mdl_rudd_day$Description[15],
                                               rdm_mdl_rudd_day$Description[16],
                                               rdm_mdl_rudd_day$Description[17],
                                               rdm_mdl_rudd_day$Description[18]),
                               component = c(rdm_mdl_rudd_day$Component[1],
                                             rdm_mdl_rudd_day$Component[2],
                                             rdm_mdl_rudd_day$Component[3],
                                             rdm_mdl_rudd_day$Component[4],
                                             rdm_mdl_rudd_day$Component[5],
                                             rdm_mdl_rudd_day$Component[6],
                                             rdm_mdl_rudd_day$Component[7],
                                             rdm_mdl_rudd_day$Component[8],
                                             rdm_mdl_rudd_day$Component[9],
                                             rdm_mdl_rudd_day$Component[10],
                                             rdm_mdl_rudd_day$Component[11],
                                             rdm_mdl_rudd_day$Component[12],
                                             rdm_mdl_rudd_day$Component[13],
                                             rdm_mdl_rudd_day$Component[14],
                                             rdm_mdl_rudd_day$Component[15],
                                             rdm_mdl_rudd_day$Component[16],
                                             rdm_mdl_rudd_day$Component[17],
                                             rdm_mdl_rudd_day$Component[18]),
                               type = c(rdm_mdl_rudd_day$Type[1],
                                        rdm_mdl_rudd_day$Type[2],
                                        rdm_mdl_rudd_day$Type[3],
                                        rdm_mdl_rudd_day$Type[4],
                                        rdm_mdl_rudd_day$Type[5],
                                        rdm_mdl_rudd_day$Type[6],
                                        rdm_mdl_rudd_day$Type[7],
                                        rdm_mdl_rudd_day$Type[8],
                                        rdm_mdl_rudd_day$Type[9],
                                        rdm_mdl_rudd_day$Type[10],
                                        rdm_mdl_rudd_day$Type[11],
                                        rdm_mdl_rudd_day$Type[12],
                                        rdm_mdl_rudd_day$Type[13],
                                        rdm_mdl_rudd_day$Type[14],
                                        rdm_mdl_rudd_day$Type[15],
                                        rdm_mdl_rudd_day$Type[16],
                                        rdm_mdl_rudd_day$Type[17],
                                        rdm_mdl_rudd_day$Type[18]),
                               term = c(rdm_mdl_rudd_day$Term[1],
                                        rdm_mdl_rudd_day$Term[2],
                                        rdm_mdl_rudd_day$Term[3],
                                        rdm_mdl_rudd_day$Term[4],
                                        rdm_mdl_rudd_day$Term[5],
                                        rdm_mdl_rudd_day$Term[6],
                                        rdm_mdl_rudd_day$Term[7],
                                        rdm_mdl_rudd_day$Term[8],
                                        rdm_mdl_rudd_day$Term[9],
                                        rdm_mdl_rudd_day$Term[10],
                                        rdm_mdl_rudd_day$Term[11],
                                        rdm_mdl_rudd_day$Term[12],
                                        rdm_mdl_rudd_day$Term[13],
                                        rdm_mdl_rudd_day$Term[14],
                                        rdm_mdl_rudd_day$Term[15],
                                        rdm_mdl_rudd_day$Term[16],
                                        rdm_mdl_rudd_day$Term[17],
                                        rdm_mdl_rudd_day$Term[18]),
                               value = c(rdm_mdl_rudd_day$Value[1],
                                         rdm_mdl_rudd_day$Value[2],
                                         rdm_mdl_rudd_day$Value[3],
                                         rdm_mdl_rudd_day$Value[4],
                                         rdm_mdl_rudd_day$Value[5],
                                         rdm_mdl_rudd_day$Value[6],
                                         rdm_mdl_rudd_day$Value[7],
                                         rdm_mdl_rudd_day$Value[8],
                                         rdm_mdl_rudd_day$Value[9],
                                         rdm_mdl_rudd_day$Value[10],
                                         rdm_mdl_rudd_day$Value[11],
                                         rdm_mdl_rudd_day$Value[12],
                                         rdm_mdl_rudd_day$Value[13],
                                         rdm_mdl_rudd_day$Value[14],
                                         rdm_mdl_rudd_day$Value[15],
                                         rdm_mdl_rudd_day$Value[16],
                                         rdm_mdl_rudd_day$Value[17],
                                         rdm_mdl_rudd_day$Value[18]))

#Percentage explained by each variable - Rudd - Day ####
#Sum of squares to give the percentage explained by each variable in the model
#Day
mdl_random_slo_int_rudd_day_aov <- anova(mdl_random_slo_int_day_rudd)
rudd_mdl_day_percentage <- tibble(species = "rudd",
                                   diel_period = "day",
                                   variable = c('Thermocline thickness (m)',
                                                'Thermocline strength (ºC)',
                                                'Thermocline depth (m)',
                                                'Seiche strength'),
                                   percentage_explained = c(
                                     round(mdl_random_slo_int_rudd_day_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_rudd_day_aov$`Sum Sq`),2),
                                     round(mdl_random_slo_int_rudd_day_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_rudd_day_aov$`Sum Sq`),2),
                                     round(mdl_random_slo_int_rudd_day_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_rudd_day_aov$`Sum Sq`),2),
                                     round(mdl_random_slo_int_rudd_day_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_rudd_day_aov$`Sum Sq`),2)
                                   )
)

rudd_day_mdl_summary <- list(fixed_effects = mdl_rudd_day_fixed,
                              percentage_explained = rudd_mdl_day_percentage,
                              random_effects = mdl_rudd_day_random)

#Night
mdl_random_slo_int_rudd_night_aov <-  anova(mdl_random_slo_int_night_rudd)
rudd_mdl_night_percentage <- tibble(species = "rudd",
                                     diel_period = "night",
                                     variable = c('Thermocline thickness (m)',
                                                  'Thermocline strength (ºC)',
                                                  'Thermocline depth (m)',
                                                  'Seiche strength'),
                                     percentage_explained = c(
                                       round(mdl_random_slo_int_rudd_night_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_rudd_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_rudd_night_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_rudd_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_rudd_night_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_rudd_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_rudd_night_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_rudd_night_aov$`Sum Sq`),2))
)

#Rudd (day & night)
rudd_mld_percentages <- bind_rows(rudd_mdl_day_percentage, rudd_mdl_night_percentage)

#Fixed effects summary - Rudd - Night ####
params_mld_rudd_night <- parameters::model_parameters(mdl_random_slo_int_night_rudd)

mdl_rudd_night_fixed <- tibble(species = "rudd", 
                                diel_period = "night",
                                parameter = c(params_mld_rudd_night$Parameter[1],
                                              params_mld_rudd_night$Parameter[2],
                                              params_mld_rudd_night$Parameter[3],
                                              params_mld_rudd_night$Parameter[4],
                                              params_mld_rudd_night$Parameter[5]),
                                coefficient = c(params_mld_rudd_night$Coefficient[1],
                                                params_mld_rudd_night$Coefficient[2],
                                                params_mld_rudd_night$Coefficient[3],
                                                params_mld_rudd_night$Coefficient[4],
                                                params_mld_rudd_night$Coefficient[5]),
                                se = c(params_mld_rudd_night$SE[1],
                                       params_mld_rudd_night$SE[2],
                                       params_mld_rudd_night$SE[3],
                                       params_mld_rudd_night$SE[4],
                                       params_mld_rudd_night$SE[5]),
                                ci_05 = c(params_mld_rudd_night$CI_low[1],
                                          params_mld_rudd_night$CI_low[2],
                                          params_mld_rudd_night$CI_low[3],
                                          params_mld_rudd_night$CI_low[4],
                                          params_mld_rudd_night$CI_low[5]),
                                ci_95 = c(params_mld_rudd_night$CI_high[1],
                                          params_mld_rudd_night$CI_high[2],
                                          params_mld_rudd_night$CI_high[3],
                                          params_mld_rudd_night$CI_high[4],
                                          params_mld_rudd_night$CI_high[5]),
                                t_value = c(params_mld_rudd_night$t[1],
                                            params_mld_rudd_night$t[2],
                                            params_mld_rudd_night$t[3],
                                            params_mld_rudd_night$t[4],
                                            params_mld_rudd_night$t[5]),
                                p_value = c(params_mld_rudd_night$p[1],
                                            params_mld_rudd_night$p[2],
                                            params_mld_rudd_night$p[3],
                                            params_mld_rudd_night$p[4],
                                            params_mld_rudd_night$p[5]))

#Random effects summary - Rudd - Night ####
rdm_mdl_rudd_night <- parameters::random_parameters(mdl_random_slo_int_night_rudd)

mdl_rudd_night_random <- tibble(species = "rudd",
                                 diel_period = "night",
                                 description = c(rdm_mdl_rudd_night$Description[1],
                                                 rdm_mdl_rudd_night$Description[2],
                                                 rdm_mdl_rudd_night$Description[3],
                                                 rdm_mdl_rudd_night$Description[4],
                                                 rdm_mdl_rudd_night$Description[5],
                                                 rdm_mdl_rudd_night$Description[6],
                                                 rdm_mdl_rudd_night$Description[7],
                                                 rdm_mdl_rudd_night$Description[8],
                                                 rdm_mdl_rudd_night$Description[9],
                                                 rdm_mdl_rudd_night$Description[10],
                                                 rdm_mdl_rudd_night$Description[11],
                                                 rdm_mdl_rudd_night$Description[12],
                                                 rdm_mdl_rudd_night$Description[13],
                                                 rdm_mdl_rudd_night$Description[14],
                                                 rdm_mdl_rudd_night$Description[15],
                                                 rdm_mdl_rudd_night$Description[16],
                                                 rdm_mdl_rudd_night$Description[17],
                                                 rdm_mdl_rudd_night$Description[18]),
                                 component = c(rdm_mdl_rudd_night$Component[1],
                                               rdm_mdl_rudd_night$Component[2],
                                               rdm_mdl_rudd_night$Component[3],
                                               rdm_mdl_rudd_night$Component[4],
                                               rdm_mdl_rudd_night$Component[5],
                                               rdm_mdl_rudd_night$Component[6],
                                               rdm_mdl_rudd_night$Component[7],
                                               rdm_mdl_rudd_night$Component[8],
                                               rdm_mdl_rudd_night$Component[9],
                                               rdm_mdl_rudd_night$Component[10],
                                               rdm_mdl_rudd_night$Component[11],
                                               rdm_mdl_rudd_night$Component[12],
                                               rdm_mdl_rudd_night$Component[13],
                                               rdm_mdl_rudd_night$Component[14],
                                               rdm_mdl_rudd_night$Component[15],
                                               rdm_mdl_rudd_night$Component[16],
                                               rdm_mdl_rudd_night$Component[17],
                                               rdm_mdl_rudd_night$Component[18]),
                                 type = c(rdm_mdl_rudd_night$Type[1],
                                          rdm_mdl_rudd_night$Type[2],
                                          rdm_mdl_rudd_night$Type[3],
                                          rdm_mdl_rudd_night$Type[4],
                                          rdm_mdl_rudd_night$Type[5],
                                          rdm_mdl_rudd_night$Type[6],
                                          rdm_mdl_rudd_night$Type[7],
                                          rdm_mdl_rudd_night$Type[8],
                                          rdm_mdl_rudd_night$Type[9],
                                          rdm_mdl_rudd_night$Type[10],
                                          rdm_mdl_rudd_night$Type[11],
                                          rdm_mdl_rudd_night$Type[12],
                                          rdm_mdl_rudd_night$Type[13],
                                          rdm_mdl_rudd_night$Type[14],
                                          rdm_mdl_rudd_night$Type[15],
                                          rdm_mdl_rudd_night$Type[16],
                                          rdm_mdl_rudd_night$Type[17],
                                          rdm_mdl_rudd_night$Type[18]),
                                 term = c(rdm_mdl_rudd_night$Term[1],
                                          rdm_mdl_rudd_night$Term[2],
                                          rdm_mdl_rudd_night$Term[3],
                                          rdm_mdl_rudd_night$Term[4],
                                          rdm_mdl_rudd_night$Term[5],
                                          rdm_mdl_rudd_night$Term[6],
                                          rdm_mdl_rudd_night$Term[7],
                                          rdm_mdl_rudd_night$Term[8],
                                          rdm_mdl_rudd_night$Term[9],
                                          rdm_mdl_rudd_night$Term[10],
                                          rdm_mdl_rudd_night$Term[11],
                                          rdm_mdl_rudd_night$Term[12],
                                          rdm_mdl_rudd_night$Term[13],
                                          rdm_mdl_rudd_night$Term[14],
                                          rdm_mdl_rudd_night$Term[15],
                                          rdm_mdl_rudd_night$Term[16],
                                          rdm_mdl_rudd_night$Term[17],
                                          rdm_mdl_rudd_night$Term[18]),
                                 value = c(rdm_mdl_rudd_night$Value[1],
                                           rdm_mdl_rudd_night$Value[2],
                                           rdm_mdl_rudd_night$Value[3],
                                           rdm_mdl_rudd_night$Value[4],
                                           rdm_mdl_rudd_night$Value[5],
                                           rdm_mdl_rudd_night$Value[6],
                                           rdm_mdl_rudd_night$Value[7],
                                           rdm_mdl_rudd_night$Value[8],
                                           rdm_mdl_rudd_night$Value[9],
                                           rdm_mdl_rudd_night$Value[10],
                                           rdm_mdl_rudd_night$Value[11],
                                           rdm_mdl_rudd_night$Value[12],
                                           rdm_mdl_rudd_night$Value[13],
                                           rdm_mdl_rudd_night$Value[14],
                                           rdm_mdl_rudd_night$Value[15],
                                           rdm_mdl_rudd_night$Value[16],
                                           rdm_mdl_rudd_night$Value[17],
                                           rdm_mdl_rudd_night$Value[18]))
#Percentage explained by each variable - Rudd - Night ####
#Night
rudd_mdl_night_percentage <- tibble(species = "rudd",
                                     diel_period = "night",
                                     variable = c('Thermocline thickness (m)',
                                                  'Thermocline strength (ºC)',
                                                  'Thermocline depth (m)',
                                                  'Seiche strength'),
                                     percentage_explained = c(
                                       round(mdl_random_slo_int_rudd_night_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_rudd_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_rudd_night_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_rudd_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_rudd_night_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_rudd_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_rudd_night_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_rudd_night_aov$`Sum Sq`),2))
)


rudd_night_mdl_summary <- list(fixed_effects = mdl_rudd_night_fixed,
                                percentage_explained = rudd_mdl_night_percentage,
                                random_effects = mdl_rudd_night_random)

#Summary of the rudd models - Rudd####
rudd_mld_summary <- list(fixed_effects = bind_rows(mdl_rudd_day_fixed,
                                                    mdl_rudd_night_fixed),
                          random_effects = bind_rows(mdl_rudd_day_random,
                                                     mdl_rudd_night_random),
                          percentage_explained = bind_rows(rudd_mdl_day_percentage,
                                                           rudd_mdl_night_percentage))

#Roach ####
#Models run - Roach ####
#random intercept and random slope per fish id - Day (Seiche only)
tic('It took this amount of time to run the model')
mdl_random_slo_int_day_roach <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + 
                                      det_therm_strength + lake_therm_depth_smoothed_center + 
                                      det_therm_deviation_center + 
                                      (1 + lake_therm_thickness_smoothed + 
                                         det_therm_strength + 
                                         lake_therm_depth_smoothed_center +
                                         det_therm_deviation_center | fishid),
                                    data = detections %>% 
                                      filter(species == "roach" & 
                                               diel_period == 'day' & 
                                               is_valid_seiche == TRUE), 
                                    REML = T, 
                                    lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()

#random intercept and random slope per fish id - Night (Seiche only)
tic('It took this amount of time to run the model')
mdl_random_slo_int_night_roach <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + 
                                        det_therm_strength + lake_therm_depth_smoothed_center + 
                                        det_therm_deviation_center + 
                                        (1 + lake_therm_thickness_smoothed + 
                                           det_therm_strength + 
                                           lake_therm_depth_smoothed_center + 
                                           det_therm_deviation_center | fishid),
                                      data = detections %>% 
                                        filter(species == "rudd" & 
                                                 diel_period == 'night' & 
                                                 is_valid_seiche == TRUE), 
                                      REML = T, 
                                      lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()

#Summary of the model results - Roach ####
tab_model(mdl_random_slo_int_day_roach, 
          mdl_random_slo_int_night_roach, dv.labels = c("Day","Night"), 
          title='Roach (Seiche only)',
          pred.labels = c('Intercept', 'Thermocline thickness (m)',
                          'Thermocline strength (ºC)',
                          'Thermocline depth (m)',
                          'Seiche strength'))#sigma² = withing group variance; tau00 = between group variance; ICC = proportion of variance explained for each group level

#Fixed effects summary - Roach - Day ####
params_mld_roach_day <- parameters::model_parameters(mdl_random_slo_int_day_roach)

mdl_roach_day_fixed <- tibble(species = "roach", 
                              diel_period = "day",
                              parameter = c(params_mld_roach_day$Parameter[1],
                                            params_mld_roach_day$Parameter[2],
                                            params_mld_roach_day$Parameter[3],
                                            params_mld_roach_day$Parameter[4],
                                            params_mld_roach_day$Parameter[5]),
                              coefficient = c(params_mld_roach_day$Coefficient[1],
                                              params_mld_roach_day$Coefficient[2],
                                              params_mld_roach_day$Coefficient[3],
                                              params_mld_roach_day$Coefficient[4],
                                              params_mld_roach_day$Coefficient[5]),
                              se = c(params_mld_roach_day$SE[1],
                                     params_mld_roach_day$SE[2],
                                     params_mld_roach_day$SE[3],
                                     params_mld_roach_day$SE[4],
                                     params_mld_roach_day$SE[5]),
                              ci_05 = c(params_mld_roach_day$CI_low[1],
                                        params_mld_roach_day$CI_low[2],
                                        params_mld_roach_day$CI_low[3],
                                        params_mld_roach_day$CI_low[4],
                                        params_mld_roach_day$CI_low[5]),
                              ci_95 = c(params_mld_roach_day$CI_high[1],
                                        params_mld_roach_day$CI_high[2],
                                        params_mld_roach_day$CI_high[3],
                                        params_mld_roach_day$CI_high[4],
                                        params_mld_roach_day$CI_high[5]),
                              t_value = c(params_mld_roach_day$t[1],
                                          params_mld_roach_day$t[2],
                                          params_mld_roach_day$t[3],
                                          params_mld_roach_day$t[4],
                                          params_mld_roach_day$t[5]),
                              p_value = c(params_mld_roach_day$p[1],
                                          params_mld_roach_day$p[2],
                                          params_mld_roach_day$p[3],
                                          params_mld_roach_day$p[4],
                                          params_mld_roach_day$p[5]))

#Random effects summary - Roach - Day ####
rdm_mdl_roach_day <- parameters::random_parameters(mdl_random_slo_int_day_roach)

mdl_roach_day_random <- tibble(species = "roach",
                               diel_period = "day",
                               description = c(rdm_mdl_roach_day$Description[1],
                                               rdm_mdl_roach_day$Description[2],
                                               rdm_mdl_roach_day$Description[3],
                                               rdm_mdl_roach_day$Description[4],
                                               rdm_mdl_roach_day$Description[5],
                                               rdm_mdl_roach_day$Description[6],
                                               rdm_mdl_roach_day$Description[7],
                                               rdm_mdl_roach_day$Description[8],
                                               rdm_mdl_roach_day$Description[9],
                                               rdm_mdl_roach_day$Description[10],
                                               rdm_mdl_roach_day$Description[11],
                                               rdm_mdl_roach_day$Description[12],
                                               rdm_mdl_roach_day$Description[13],
                                               rdm_mdl_roach_day$Description[14],
                                               rdm_mdl_roach_day$Description[15],
                                               rdm_mdl_roach_day$Description[16],
                                               rdm_mdl_roach_day$Description[17],
                                               rdm_mdl_roach_day$Description[18]),
                               component = c(rdm_mdl_roach_day$Component[1],
                                             rdm_mdl_roach_day$Component[2],
                                             rdm_mdl_roach_day$Component[3],
                                             rdm_mdl_roach_day$Component[4],
                                             rdm_mdl_roach_day$Component[5],
                                             rdm_mdl_roach_day$Component[6],
                                             rdm_mdl_roach_day$Component[7],
                                             rdm_mdl_roach_day$Component[8],
                                             rdm_mdl_roach_day$Component[9],
                                             rdm_mdl_roach_day$Component[10],
                                             rdm_mdl_roach_day$Component[11],
                                             rdm_mdl_roach_day$Component[12],
                                             rdm_mdl_roach_day$Component[13],
                                             rdm_mdl_roach_day$Component[14],
                                             rdm_mdl_roach_day$Component[15],
                                             rdm_mdl_roach_day$Component[16],
                                             rdm_mdl_roach_day$Component[17],
                                             rdm_mdl_roach_day$Component[18]),
                               type = c(rdm_mdl_roach_day$Type[1],
                                        rdm_mdl_roach_day$Type[2],
                                        rdm_mdl_roach_day$Type[3],
                                        rdm_mdl_roach_day$Type[4],
                                        rdm_mdl_roach_day$Type[5],
                                        rdm_mdl_roach_day$Type[6],
                                        rdm_mdl_roach_day$Type[7],
                                        rdm_mdl_roach_day$Type[8],
                                        rdm_mdl_roach_day$Type[9],
                                        rdm_mdl_roach_day$Type[10],
                                        rdm_mdl_roach_day$Type[11],
                                        rdm_mdl_roach_day$Type[12],
                                        rdm_mdl_roach_day$Type[13],
                                        rdm_mdl_roach_day$Type[14],
                                        rdm_mdl_roach_day$Type[15],
                                        rdm_mdl_roach_day$Type[16],
                                        rdm_mdl_roach_day$Type[17],
                                        rdm_mdl_roach_day$Type[18]),
                               term = c(rdm_mdl_roach_day$Term[1],
                                        rdm_mdl_roach_day$Term[2],
                                        rdm_mdl_roach_day$Term[3],
                                        rdm_mdl_roach_day$Term[4],
                                        rdm_mdl_roach_day$Term[5],
                                        rdm_mdl_roach_day$Term[6],
                                        rdm_mdl_roach_day$Term[7],
                                        rdm_mdl_roach_day$Term[8],
                                        rdm_mdl_roach_day$Term[9],
                                        rdm_mdl_roach_day$Term[10],
                                        rdm_mdl_roach_day$Term[11],
                                        rdm_mdl_roach_day$Term[12],
                                        rdm_mdl_roach_day$Term[13],
                                        rdm_mdl_roach_day$Term[14],
                                        rdm_mdl_roach_day$Term[15],
                                        rdm_mdl_roach_day$Term[16],
                                        rdm_mdl_roach_day$Term[17],
                                        rdm_mdl_roach_day$Term[18]),
                               value = c(rdm_mdl_roach_day$Value[1],
                                         rdm_mdl_roach_day$Value[2],
                                         rdm_mdl_roach_day$Value[3],
                                         rdm_mdl_roach_day$Value[4],
                                         rdm_mdl_roach_day$Value[5],
                                         rdm_mdl_roach_day$Value[6],
                                         rdm_mdl_roach_day$Value[7],
                                         rdm_mdl_roach_day$Value[8],
                                         rdm_mdl_roach_day$Value[9],
                                         rdm_mdl_roach_day$Value[10],
                                         rdm_mdl_roach_day$Value[11],
                                         rdm_mdl_roach_day$Value[12],
                                         rdm_mdl_roach_day$Value[13],
                                         rdm_mdl_roach_day$Value[14],
                                         rdm_mdl_roach_day$Value[15],
                                         rdm_mdl_roach_day$Value[16],
                                         rdm_mdl_roach_day$Value[17],
                                         rdm_mdl_roach_day$Value[18]))

#Sum of squares to give the percentage explained by each variable in the model
#Day
mdl_random_slo_int_roach_day_aov <-  anova(mdl_random_slo_int_day_roach)
#Night
mdl_random_slo_int_roach_night_aov <-  anova(mdl_random_slo_int_night_roach)

#Percentage explained by each variable - Roach - Day ####
#Day
roach_mdl_day_percentage <- tibble(species = "roach",
                                  diel_period = "day",
                                  variable = c('Thermocline thickness (m)',
                                               'Thermocline strength (ºC)',
                                               'Thermocline depth (m)',
                                               'Seiche strength'),
                                  percentage_explained = c(
                                    round(mdl_random_slo_int_roach_day_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_roach_day_aov$`Sum Sq`),2),
                                    round(mdl_random_slo_int_roach_day_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_roach_day_aov$`Sum Sq`),2),
                                    round(mdl_random_slo_int_roach_day_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_roach_day_aov$`Sum Sq`),2),
                                    round(mdl_random_slo_int_roach_day_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_roach_day_aov$`Sum Sq`),2)
                                  )
)

roach_day_mdl_summary <- list(fixed_effects = mdl_roach_day_fixed,
                              percentage_explained = roach_mdl_day_percentage,
                              random_effects = mdl_roach_day_random)


#Fixed effects summary - Roach - Night ####
params_mld_roach_night <- parameters::model_parameters(mdl_random_slo_int_night_roach)

mdl_roach_night_fixed <- tibble(species = "roach", 
                              diel_period = "night",
                              parameter = c(params_mld_roach_night$Parameter[1],
                                            params_mld_roach_night$Parameter[2],
                                            params_mld_roach_night$Parameter[3],
                                            params_mld_roach_night$Parameter[4],
                                            params_mld_roach_night$Parameter[5]),
                              coefficient = c(params_mld_roach_night$Coefficient[1],
                                              params_mld_roach_night$Coefficient[2],
                                              params_mld_roach_night$Coefficient[3],
                                              params_mld_roach_night$Coefficient[4],
                                              params_mld_roach_night$Coefficient[5]),
                              se = c(params_mld_roach_night$SE[1],
                                     params_mld_roach_night$SE[2],
                                     params_mld_roach_night$SE[3],
                                     params_mld_roach_night$SE[4],
                                     params_mld_roach_night$SE[5]),
                              ci_05 = c(params_mld_roach_night$CI_low[1],
                                        params_mld_roach_night$CI_low[2],
                                        params_mld_roach_night$CI_low[3],
                                        params_mld_roach_night$CI_low[4],
                                        params_mld_roach_night$CI_low[5]),
                              ci_95 = c(params_mld_roach_night$CI_high[1],
                                        params_mld_roach_night$CI_high[2],
                                        params_mld_roach_night$CI_high[3],
                                        params_mld_roach_night$CI_high[4],
                                        params_mld_roach_night$CI_high[5]),
                              t_value = c(params_mld_roach_night$t[1],
                                          params_mld_roach_night$t[2],
                                          params_mld_roach_night$t[3],
                                          params_mld_roach_night$t[4],
                                          params_mld_roach_night$t[5]),
                              p_value = c(params_mld_roach_night$p[1],
                                          params_mld_roach_night$p[2],
                                          params_mld_roach_night$p[3],
                                          params_mld_roach_night$p[4],
                                          params_mld_roach_night$p[5]))

#Random effects summary - Roach - Night ####
rdm_mdl_roach_night <- parameters::random_parameters(mdl_random_slo_int_night_roach)

mdl_roach_night_random <- tibble(species = "roach",
                               diel_period = "night",
                               description = c(rdm_mdl_roach_night$Description[1],
                                               rdm_mdl_roach_night$Description[2],
                                               rdm_mdl_roach_night$Description[3],
                                               rdm_mdl_roach_night$Description[4],
                                               rdm_mdl_roach_night$Description[5],
                                               rdm_mdl_roach_night$Description[6],
                                               rdm_mdl_roach_night$Description[7],
                                               rdm_mdl_roach_night$Description[8],
                                               rdm_mdl_roach_night$Description[9],
                                               rdm_mdl_roach_night$Description[10],
                                               rdm_mdl_roach_night$Description[11],
                                               rdm_mdl_roach_night$Description[12],
                                               rdm_mdl_roach_night$Description[13],
                                               rdm_mdl_roach_night$Description[14],
                                               rdm_mdl_roach_night$Description[15],
                                               rdm_mdl_roach_night$Description[16],
                                               rdm_mdl_roach_night$Description[17],
                                               rdm_mdl_roach_night$Description[18]),
                               component = c(rdm_mdl_roach_night$Component[1],
                                             rdm_mdl_roach_night$Component[2],
                                             rdm_mdl_roach_night$Component[3],
                                             rdm_mdl_roach_night$Component[4],
                                             rdm_mdl_roach_night$Component[5],
                                             rdm_mdl_roach_night$Component[6],
                                             rdm_mdl_roach_night$Component[7],
                                             rdm_mdl_roach_night$Component[8],
                                             rdm_mdl_roach_night$Component[9],
                                             rdm_mdl_roach_night$Component[10],
                                             rdm_mdl_roach_night$Component[11],
                                             rdm_mdl_roach_night$Component[12],
                                             rdm_mdl_roach_night$Component[13],
                                             rdm_mdl_roach_night$Component[14],
                                             rdm_mdl_roach_night$Component[15],
                                             rdm_mdl_roach_night$Component[16],
                                             rdm_mdl_roach_night$Component[17],
                                             rdm_mdl_roach_night$Component[18]),
                               type = c(rdm_mdl_roach_night$Type[1],
                                        rdm_mdl_roach_night$Type[2],
                                        rdm_mdl_roach_night$Type[3],
                                        rdm_mdl_roach_night$Type[4],
                                        rdm_mdl_roach_night$Type[5],
                                        rdm_mdl_roach_night$Type[6],
                                        rdm_mdl_roach_night$Type[7],
                                        rdm_mdl_roach_night$Type[8],
                                        rdm_mdl_roach_night$Type[9],
                                        rdm_mdl_roach_night$Type[10],
                                        rdm_mdl_roach_night$Type[11],
                                        rdm_mdl_roach_night$Type[12],
                                        rdm_mdl_roach_night$Type[13],
                                        rdm_mdl_roach_night$Type[14],
                                        rdm_mdl_roach_night$Type[15],
                                        rdm_mdl_roach_night$Type[16],
                                        rdm_mdl_roach_night$Type[17],
                                        rdm_mdl_roach_night$Type[18]),
                               term = c(rdm_mdl_roach_night$Term[1],
                                        rdm_mdl_roach_night$Term[2],
                                        rdm_mdl_roach_night$Term[3],
                                        rdm_mdl_roach_night$Term[4],
                                        rdm_mdl_roach_night$Term[5],
                                        rdm_mdl_roach_night$Term[6],
                                        rdm_mdl_roach_night$Term[7],
                                        rdm_mdl_roach_night$Term[8],
                                        rdm_mdl_roach_night$Term[9],
                                        rdm_mdl_roach_night$Term[10],
                                        rdm_mdl_roach_night$Term[11],
                                        rdm_mdl_roach_night$Term[12],
                                        rdm_mdl_roach_night$Term[13],
                                        rdm_mdl_roach_night$Term[14],
                                        rdm_mdl_roach_night$Term[15],
                                        rdm_mdl_roach_night$Term[16],
                                        rdm_mdl_roach_night$Term[17],
                                        rdm_mdl_roach_night$Term[18]),
                               value = c(rdm_mdl_roach_night$Value[1],
                                         rdm_mdl_roach_night$Value[2],
                                         rdm_mdl_roach_night$Value[3],
                                         rdm_mdl_roach_night$Value[4],
                                         rdm_mdl_roach_night$Value[5],
                                         rdm_mdl_roach_night$Value[6],
                                         rdm_mdl_roach_night$Value[7],
                                         rdm_mdl_roach_night$Value[8],
                                         rdm_mdl_roach_night$Value[9],
                                         rdm_mdl_roach_night$Value[10],
                                         rdm_mdl_roach_night$Value[11],
                                         rdm_mdl_roach_night$Value[12],
                                         rdm_mdl_roach_night$Value[13],
                                         rdm_mdl_roach_night$Value[14],
                                         rdm_mdl_roach_night$Value[15],
                                         rdm_mdl_roach_night$Value[16],
                                         rdm_mdl_roach_night$Value[17],
                                         rdm_mdl_roach_night$Value[18]))

#Percentage explained by each variable - Roach - Night ####
#Night
roach_mdl_night_percentage <- tibble(species = "roach",
                                     diel_period = "night",
                                     variable = c('Thermocline thickness (m)',
                                                  'Thermocline strength (ºC)',
                                                  'Thermocline depth (m)',
                                                  'Seiche strength'),
                                     percentage_explained = c(
                                       round(mdl_random_slo_int_roach_night_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_roach_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_roach_night_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_roach_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_roach_night_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_roach_night_aov$`Sum Sq`),2),
                                       round(mdl_random_slo_int_roach_night_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_roach_night_aov$`Sum Sq`),2))
)


roach_night_mdl_summary <- list(fixed_effects = mdl_roach_night_fixed,
                              percentage_explained = roach_mdl_night_percentage,
                              random_effects = mdl_roach_night_random)

#Summary of the roach models - Roach ####
roach_mld_summary <- list(fixed_effects = bind_rows(mdl_roach_day_fixed,
                                                    mdl_roach_night_fixed),
                          random_effects = bind_rows(mdl_roach_day_random,
                                                     mdl_roach_night_random),
                          percentage_explained = bind_rows(roach_mdl_day_percentage,
                                                           roach_mdl_night_percentage))

#Roach (day & night)
roach_mld_percentages <- bind_rows(roach_mdl_day_percentage, roach_mdl_night_percentage)

mld_percentages <- bind_rows(pike_mld_percentages, 
                             wels_mld_percentages,
                             tench_mld_percentages,
                             rudd_mld_percentages,
                             roach_mld_percentages)

#Binding all models results in one object####
all_mlds_results <- list(fixed_effects = bind_rows(pike_mld_summary$fixed_effects,
                                 wels_mld_summary$fixed_effects,
                                 tench_mld_summary$fixed_effects,
                                 rudd_mld_summary$fixed_effects,
                                 roach_mld_summary$fixed_effects),
       random_effects = bind_rows(pike_mld_summary$random_effects,
                                  wels_mld_summary$random_effects,
                                  tench_mld_summary$random_effects,
                                  rudd_mld_summary$random_effects,
                                  roach_mld_summary$random_effects),
       percentage_explained = bind_rows(pike_mld_summary$percentage_explained,
                                        wels_mld_summary$percentage_explained,
                                        tench_mld_summary$percentage_explained,
                                        rudd_mld_summary$percentage_explained,
                                        roach_mld_summary$percentage_explained))

#Saving the model results
#Fixed effects
xlsx::write.xlsx(x = all_mlds_results$fixed_effects, 
                 sheetName = "Fixed effects", 
                 file = here::here('data', 'products', 'models_summary.xlsx'),
                 append = F)
#Random effects
xlsx::write.xlsx(x = all_mlds_results$random_effects, 
                 sheetName = "Random effects", 
                 file = here::here('data', 'products', 'models_summary.xlsx'),
                 append = T)
#Percentage explained
xlsx::write.xlsx(x = all_mlds_results$percentage_explained, 
                 sheetName = "Percentage explained", 
                 file = here::here('data', 'products', 'models_summary.xlsx'),
                 append = T)

# write_csv(x = mdl_percentages, 
#           file = here::here('data', 'products', 'mdl_variables_effects_percentages.csv'))

#Quick dataviz
mld_percentages %>%
  mutate(diel_period = fct_recode(.f = diel_period, 
                                  "Day" = "day", "Night" = "night")) %>%
  mutate(species = fct_recode(.f = species, 
                              "Pike" = "pike", "Roach" = "roach",
                              "Rudd" = "rudd", "Tench" = "tench", 
                              "Wels catfish" = "wels")) %>%
  ggplot(aes(x = diel_period, y = percentage_explained, 
             color = variable, fill = variable)) + 
  facet_wrap(~ species) +
  geom_bar(position="fill", stat="identity")+
  labs(title = "Linear Mixed Effects model", subtitle = "Milada lake",
       x = "Diel period", y = "% explained", fill = "", color = "") +
  scale_color_viridis_d(option = "C")+
  scale_fill_viridis_d(option = "C")+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggsave(filename = here::here('outputs', 'models_percent_explained.jpeg'), 
         device = "jpeg", units = "cm", dpi = "retina", width = 18, height = 12)


# ##Legacy models ####
# #lme4 package version
# #global_mdl2 <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + diel_period + EaWe_diff + lake_therm_depth_smoothed_center + (1 + lake_therm_thickness_smoothed + diel_period + EaWe_diff + lake_therm_depth_smoothed_center|fishid),
# #                    data = detections_mdl, REML = F) #intercept varying among years and among fish_id within years (nested random effects)
# 
# #lme4 syntax
# #Random intercept per fishid
# tic("model fitting")
# global_mdl1 <- lmer(formula = det_depth ~ diel_period + lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1|fishid),
#                     data = detections_mdl, REML = T, 
#                     lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
# toc()
# 
# # tic("model fitting")
# # global_mdl1_thick <- update(global_mdl1, . ~ . - lake_therm_thickness_smoothed)
# # toc()
# 
# ##Random intercept per fishid - Day####
# tic("model fitting")
# global_mdl_inter_day <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1|fishid),
#                              data = detections_mdl[diel_period=='day'], REML = T, 
#                              lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
# toc()
# 
# ##Random intercept per fishid - Night####
# tic("model fitting")
# global_mdl_inter_night <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1|fishid),
#                                data = detections_mdl[diel_period=='night'], REML = T, 
#                                lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
# toc()
# 
# 
# #Table comparing day and night####
# tab_model(global_mdl_inter_day, global_mdl_inter_night, dv.labels = c("Day","Night"), 
#           title='Rudd', p.style='asterisk',
#           pred.labels = c('Intercept', 'Thermocline thickness (m)',
#                           'Thermocline strength (ºC)',
#                           'Thermocline depth (m)',
#                           'Seiche strength'))#sigma² = withing group variance; tau00 = between group variance; ICC = proportion of variance explained for each group level
# 
# #Sum of squares to give the percentage explained by each variable in the model####
# #Day
# global_mdl_inter_day_aov <-  anova(global_mdl_inter_day)
# 
# #Percentage explained by each variable
# round(global_mdl_inter_day_aov$`Sum Sq`[1]*100/sum(global_mdl_inter_day_aov$`Sum Sq`),2)
# round(global_mdl_inter_day_aov$`Sum Sq`[2]*100/sum(global_mdl_inter_day_aov$`Sum Sq`),2)
# round(global_mdl_inter_day_aov$`Sum Sq`[3]*100/sum(global_mdl_inter_day_aov$`Sum Sq`),2)
# round(global_mdl_inter_day_aov$`Sum Sq`[4]*100/sum(global_mdl_inter_day_aov$`Sum Sq`),2)
# 
# #Night
# global_mdl_inter_night_aov <-  anova(global_mdl_inter_night)
# 
# #Percentage explained by each variable
# round(global_mdl_inter_night_aov$`Sum Sq`[1]*100/sum(global_mdl_inter_night_aov$`Sum Sq`),2)
# round(global_mdl_inter_night_aov$`Sum Sq`[2]*100/sum(global_mdl_inter_night_aov$`Sum Sq`),2)
# round(global_mdl_inter_night_aov$`Sum Sq`[3]*100/sum(global_mdl_inter_night_aov$`Sum Sq`),2)
# round(global_mdl_inter_night_aov$`Sum Sq`[4]*100/sum(global_mdl_inter_night_aov$`Sum Sq`),2)
# 
# 
# ##Random intercept per fishid - Day (Seiche only)
# tic("model fitting")
# global_mdl_inter_day_seiche <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1|fishid),
#                                     data = detections_mdl[diel_period=='day' & Seiche=='on'], REML = T, 
#                                     lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
# toc()
# 
# ##Random intercept per fishid - Night (Seiche only)
# tic("model fitting")
# global_mdl_inter_night_seiche <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1|fishid),
#                                       data = detections_mdl[diel_period=='night' & Seiche=='on'], REML = T, 
#                                       lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
# toc()
# 
# #Table comparing day and night (Seiche only)####
# tab_model(global_mdl_inter_day_seiche, global_mdl_inter_night_seiche, dv.labels = c("Day","Night"), 
#           title='Rudd (Seiche only)', p.style='asterisk',
#           pred.labels = c('Intercept', 'Thermocline thickness (m)',
#                           'Thermocline strength (ºC)',
#                           'Thermocline depth (m)',
#                           'Seiche strength'))#sigma² = withing group variance; tau00 = between group variance; ICC = proportion of variance explained for each group level
# 
# 
# #Sum of squares to give the percentage explained by each variable in the model####
# #Day
# global_mdl_inter_day_Seiche_aov <-  anova(global_mdl_inter_day_seiche)
# 
# #Percentage explained by each variable
# round(global_mdl_inter_day_Seiche_aov$`Sum Sq`[1]*100/sum(global_mdl_inter_day_Seiche_aov$`Sum Sq`),2)
# round(global_mdl_inter_day_Seiche_aov$`Sum Sq`[2]*100/sum(global_mdl_inter_day_Seiche_aov$`Sum Sq`),2)
# round(global_mdl_inter_day_Seiche_aov$`Sum Sq`[3]*100/sum(global_mdl_inter_day_Seiche_aov$`Sum Sq`),2)
# round(global_mdl_inter_day_Seiche_aov$`Sum Sq`[4]*100/sum(global_mdl_inter_day_Seiche_aov$`Sum Sq`),2)
# 
# #Night
# global_mdl_inter_night_seiche_aov <-  anova(global_mdl_inter_night_seiche)
# 
# #Percentage explained by each variable
# round(global_mdl_inter_night_seiche_aov$`Sum Sq`[1]*100/sum(global_mdl_inter_night_seiche_aov$`Sum Sq`),2)
# round(global_mdl_inter_night_seiche_aov$`Sum Sq`[2]*100/sum(global_mdl_inter_night_seiche_aov$`Sum Sq`),2)
# round(global_mdl_inter_night_seiche_aov$`Sum Sq`[3]*100/sum(global_mdl_inter_night_seiche_aov$`Sum Sq`),2)
# round(global_mdl_inter_night_seiche_aov$`Sum Sq`[4]*100/sum(global_mdl_inter_night_seiche_aov$`Sum Sq`),2)
# 
# 
# 
# 
# 
# #Random intercept per fishid with interactions
# tic("model fitting")
# global_mdl2 <- lmer(formula = det_depth ~ diel_period * (lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center) + (1|fishid),
#                     data = detections_mdl, REML = T, 
#                     lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
# toc()
# 
# 
# 
# AICctab(global_mdl1, global_mdl1_thick, global_mdl2)
# 
# #Random slope per fishid
# tic("model fitting")
# global_mdl3 <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center + (1 + lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center | fishid),
#                     data = detections_mdl, REML = T,
#                     lmerControl(optimizer = 'bobyqa')) #random group slope
# toc()
# 
# 
# 
# #random intercept and random slope per fish id - Day (Seiche only) #####
# tic('Fuck, it took this amount of time to run the model')
# mdl_random_slo_int_day <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1+lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center|fishid),
#                                data = detections_mdl[diel_period=='day' & Seiche == 'on'], REML = T, 
#                                lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
# toc()
# 
# #random intercept and random slope per fish id - Night (Seiche only) #####
# tic('Fuck, it took this amount of time to run the model')
# mdl_random_slo_int_night <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1+lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center|fishid),
#                                  data = detections_mdl[diel_period=='night' & Seiche == 'on'], REML = T, 
#                                  lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
# toc()
# 
# 
# tab_model(mdl_random_slo_int_day, mdl_random_slo_int_night, dv.labels = c("Day","Night"), 
#           title='Wels (Seiche only)', p.style='asterisk',
#           pred.labels = c('Intercept', 'Thermocline thickness (m)',
#                           'Thermocline strength (ºC)',
#                           'Thermocline depth (m)',
#                           'Seiche strength'))#sigma² = withing group variance; tau00 = between group variance; ICC = proportion of variance explained for each group level
# 
# 
# 
# 
# 
# #Sum of squares to give the percentage explained by each variable in the model####
# #Day
# mdl_random_slo_int_day_aov <-  anova(mdl_random_slo_int_day)
# 
# 
# 
# #Percentage explained by each variable
# round(mdl_random_slo_int_day_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_day_aov$`Sum Sq`),2)
# round(mdl_random_slo_int_day_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_day_aov$`Sum Sq`),2)
# round(mdl_random_slo_int_day_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_day_aov$`Sum Sq`),2)
# round(mdl_random_slo_int_day_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_day_aov$`Sum Sq`),2)
# 
# #Night
# mdl_random_slo_int_night_aov <-  anova(mdl_random_slo_int_night)
# 
# #Percentage explained by each variable
# round(mdl_random_slo_int_night_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_night_aov$`Sum Sq`),2)
# round(mdl_random_slo_int_night_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_night_aov$`Sum Sq`),2)
# round(mdl_random_slo_int_night_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_night_aov$`Sum Sq`),2)
# round(mdl_random_slo_int_night_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_night_aov$`Sum Sq`),2)
# 
# plot_model(mdl_random_slo_int_day, title = 'Wels (day)')
# plot_model(mdl_random_slo_int_night, title = 'Wels (night)')
# 
# plot_model(mdl_random_slo_int_night, type = "pred", axis.lim = c(1,5),
#            terms = c('lake_therm_depth_smoothed_center', 'fishid [T412073, T412074, T412075, T412076]'))
# 
# plot_residuals(mdl_random_slo_int_night)
# 
# 
# 
# 
# #Running more than 1 model at once at different cores
# models <- c('lmer(formula = det_depth ~ lake_therm_thickness_smoothed + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center + (1|fishid), data = i, REML = F)',
#             'lmer(formula = det_depth ~ lake_therm_thickness_smoothed + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center + (1 + lake_therm_thickness_smoothed + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center | fishid), data = i, REML = F)')
# m.list <- f_lmer_mc(data = detections_mdl, calls = models, mc.cores = 20)
# 
# 
# #Diagnostic plot
# sjp.lmer(global_mdl8, type = 'fe')
# sjp.lmer(global_mdl_inter_night_seiche_aov, type = "fe.ri", sort='fishid')
# 
# plot_model(global_mdl_inter_night_seiche, type = 'est', title = 'Tench', vline.color = "black")
# 
# #Model assumptions checking####
# #homoscedasticity test
# car::leveneTest(residuals(global_mdl6) ~ detections_mdl$lake_therm_thickness_smoothed + detections_mdl$lake_therm_depth_smoothed_center + detections_mdl$diel_period + detections_mdl$det_therm_deviation_center)
# 
# #Normality of residuals
# qqnorm(residuals(global_mdl6))#OK
# 
# #Linearity in each variable
# #Check for each one of the independent variables
# ##Thermocline thickness
# ggplot(data.frame(x1= detections_mdl$lake_therm_thickness_smoothed,
#                   pearson = residuals(global_mdl6,type="pearson")), aes(x = x1, y = pearson))+
#   geom_point(alpha=0.5, shape='.')+
#   geom_smooth(method='lm', se=F)+
#   theme_bw()
# ##Balanced thermocline depth (center)
# ggplot(data.frame(x1= detections_mdl$lake_therm_depth_smoothed_center,
#                   pearson = residuals(global_mdl6,type="pearson")), aes(x = x1, y = pearson))+
#   geom_point(alpha=0.5, shape='.')+
#   geom_smooth(method='lm', se=F)+
#   theme_bw()
# ##Period of the day
# ggplot(data.frame(x1= detections_mdl$diel_period,
#                   pearson = residuals(global_mdl6,type="pearson")), aes(x = x1, y = pearson))+
#   geom_point(alpha=0.5, shape='.')+
#   geom_smooth(method='lm', se=F)+
#   theme_bw()
# ##Seiche strength
# ggplot(data.frame(x1= detections_mdl$det_therm_deviation_center,
#                   pearson = residuals(global_mdl6,type="pearson")), aes(x = x1, y = pearson))+
#   geom_point(alpha=0.5, shape='.')+
#   geom_smooth(method='lm', se=F)+
#   theme_bw()
# 
# 
# #Selection of the best model
# AICtab(global_mdl1, global_mdl2, global_mdl3, global_mdl4)#Akaike Information Criterion
# BICtab(global_mdl1, global_mdl2, global_mdl3, global_mdl4)#Bayesian Information Criterion
# 
# 
# #Inter Class Correlation (ICC) - (Measure of how much of the variation in the response variable is accounted by the random effect)
# r1Var <- as.numeric(VarCorr(global_mdl1)[["fishid"]])
# residVar <- attr(VarCorr(global_mdl1), "sc")^2
# round((r1Var/(r1Var+residVar)), 2) #percentage of variation explained by the random effect
# 
# AICtab(global_mdl8, global_mdl6)
# 
# #Model diagnosis
# ranef(global_mdl1)#estimated deviation between each level of the random effect and the overall average:
# fixef(global_mdl1)#fixed effects outputs
# round(MuMIn::r.squaredGLMM(global_mdl1),2)#R2m = R2 for the marginal test (fixed); R2c = R2 for the conditional test (fixed + random)
# 
# # #Model simplification by backwards
# # lmerTest::step(model=global_mdl2)#the analysis suggested to keep all variables into the explanatory side of the model 
# # #test the significance of each variable (full model x reduced model)
# # #Temperature
# # asp_mdl_int_temp <- lmer(formula = eggs_sca ~ dayoy_sca + visit_sca + dayNight + (1|year),
# #                          data = df_complete, REML = F) #intercept varying among years and among fish_id within years (nested random effects)
# # anova(asp_mdl_int,asp_mdl_int_temp)#Temperature is significant
# # #Period of the day
# # asp_mdl_int_day <- lmer(formula = eggs_sca ~ dayoy_sca + visit_sca + temp_sca + (1|year),
# #                         data = df_complete, REML = F) #intercept varying among years and among fish_id within years (nested random effects)
# # anova(asp_mdl_int,asp_mdl_int_day)#Period of the day is significant
# # #Bleak presence
# # asp_mdl_int_bleak <- lmer(formula = eggs_sca ~ dayoy_sca + dayNight + temp_sca + (1|year),
# #                           data = df_complete, REML = F) #intercept varying among years and among fish_id within years (nested random effects)
# # anova(asp_mdl_int,asp_mdl_int_bleak)#Bleak presence is significant
# # #day of the year
# # asp_mdl_int_dayoy <- lmer(formula = eggs_sca ~ visit_sca + dayNight + temp_sca + (1|year),
# #                           data = df_complete, REML = F) #intercept varying among years and among fish_id within years (nested random effects)
# # 
# # anova(asp_mdl_int, asp_mdl_int_dayoy)#Day of the year is significant
# 
# #stepwise selection
# MASS::stepAIC(object = global_mdl6, direction = 'backward')
# 
# #Models with each individual in separate####
# tench_ids <- unique(detections_mdl[,fishid])#vector of unique fish ids
# 
# mdl_list <- list()#creating a list to store the models
# 
# #for cycle for the computation of all models 
# for (i  in 1:length(unique(detections_mdl[,fishid]))) {
#   mdl_list[[i]] <- lm(formula = det_depth ~ therm_temp_strength + lake_therm_thickness_smoothed + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center, data = detections_mdl[fishid==tench_ids[[i]]])
# }
# 
# #apllying the summary of all models at once
# lapply(mdl_list , summary)
# 
# #creating an automatic export of the results into a table
# sjt.lm(mdl_list[[1]],mdl_list[[2]],mdl_list[[3]],mdl_list[[4]],mdl_list[[5]],
#        mdl_list[[6]],mdl_list[[7]],mdl_list[[8]],mdl_list[[9]],mdl_list[[10]],
#        mdl_list[[11]],mdl_list[[12]],mdl_list[[13]])