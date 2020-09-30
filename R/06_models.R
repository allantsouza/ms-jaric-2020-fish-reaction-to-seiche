# TODO: this should be for loop to compute for all species
species <- "pike"


# Load data ---------------------------------------------------------------
fish_raw <- read_csv(file = "data/raw/fishIDs.csv", col_types = "ccdc") %>%
  filter(species == !!species) %>%
  mutate(data_path = here("data/products/fish",paste0(tag_sn, ".csv")))


detections <- fish_raw %>%
  filter(file.exists(data_path)) %>%
  pull(data_path) %>%
  # read in all the files, appending the path before the filename
  map(~ read_csv(.)) %>% 
  reduce(rbind) %>%
  inner_join(fish_raw[,c("tag_sn","fishid")])



# Models ------------------------------------------------------------------

#Testing the data
library(lme4)
library(lmerTest)
#install.packages('MuMIn', 'v1.40.4')#older version because the newer verion is not available in R 3.5.1
library(MuMIn)
#library(devtools)
#install_version("sjPlot", version = "2.4.0", repos = "http://cran.us.r-project.org")
#install.packages('sjPlot', version = '2.6.1')
library(sjPlot)
library(bbmle)
library(hier.part)
library(parallel)
library(tictoc)#to check time elapsed
#install.packages("optimx")
library(optimx)#optmization for lmer

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




#lme4 package version
#global_mdl2 <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + diel_period + EaWe_diff + lake_therm_depth_smoothed_center + (1 + lake_therm_thickness_smoothed + diel_period + EaWe_diff + lake_therm_depth_smoothed_center|fishid),
#                    data = detections_mdl, REML = F) #intercept varying among years and among fish_id within years (nested random effects)

#lme4 syntax
#Random intercept per fishid
tic("model fitting")
global_mdl1 <- lmer(formula = det_depth ~ diel_period + lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1|fishid),
                    data = detections_mdl, REML = T, 
                    lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
toc()

# tic("model fitting")
# global_mdl1_thick <- update(global_mdl1, . ~ . - lake_therm_thickness_smoothed)
# toc()

##Random intercept per fishid - Day####
tic("model fitting")
global_mdl_inter_day <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1|fishid),
                             data = detections_mdl[diel_period=='day'], REML = T, 
                             lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
toc()

##Random intercept per fishid - Night####
tic("model fitting")
global_mdl_inter_night <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1|fishid),
                               data = detections_mdl[diel_period=='night'], REML = T, 
                               lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
toc()


#Table comparing day and night####
tab_model(global_mdl_inter_day, global_mdl_inter_night, dv.labels = c("Day","Night"), 
          title='Rudd', p.style='asterisk',
          pred.labels = c('Intercept', 'Thermocline thickness (m)',
                          'Thermocline strength (ºC)',
                          'Thermocline depth (m)',
                          'Seiche strength'))#sigma² = withing group variance; tau00 = between group variance; ICC = proportion of variance explained for each group level

#Sum of squares to give the percentage explained by each variable in the model####
#Day
global_mdl_inter_day_aov <-  anova(global_mdl_inter_day)

#Percentage explained by each variable
round(global_mdl_inter_day_aov$`Sum Sq`[1]*100/sum(global_mdl_inter_day_aov$`Sum Sq`),2)
round(global_mdl_inter_day_aov$`Sum Sq`[2]*100/sum(global_mdl_inter_day_aov$`Sum Sq`),2)
round(global_mdl_inter_day_aov$`Sum Sq`[3]*100/sum(global_mdl_inter_day_aov$`Sum Sq`),2)
round(global_mdl_inter_day_aov$`Sum Sq`[4]*100/sum(global_mdl_inter_day_aov$`Sum Sq`),2)

#Night
global_mdl_inter_night_aov <-  anova(global_mdl_inter_night)

#Percentage explained by each variable
round(global_mdl_inter_night_aov$`Sum Sq`[1]*100/sum(global_mdl_inter_night_aov$`Sum Sq`),2)
round(global_mdl_inter_night_aov$`Sum Sq`[2]*100/sum(global_mdl_inter_night_aov$`Sum Sq`),2)
round(global_mdl_inter_night_aov$`Sum Sq`[3]*100/sum(global_mdl_inter_night_aov$`Sum Sq`),2)
round(global_mdl_inter_night_aov$`Sum Sq`[4]*100/sum(global_mdl_inter_night_aov$`Sum Sq`),2)


##Random intercept per fishid - Day (Seiche only)
tic("model fitting")
global_mdl_inter_day_seiche <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1|fishid),
                                    data = detections_mdl[diel_period=='day' & Seiche=='on'], REML = T, 
                                    lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
toc()

##Random intercept per fishid - Night (Seiche only)
tic("model fitting")
global_mdl_inter_night_seiche <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1|fishid),
                                      data = detections_mdl[diel_period=='night' & Seiche=='on'], REML = T, 
                                      lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
toc()

#Table comparing day and night (Seiche only)####
tab_model(global_mdl_inter_day_seiche, global_mdl_inter_night_seiche, dv.labels = c("Day","Night"), 
          title='Rudd (Seiche only)', p.style='asterisk',
          pred.labels = c('Intercept', 'Thermocline thickness (m)',
                          'Thermocline strength (ºC)',
                          'Thermocline depth (m)',
                          'Seiche strength'))#sigma² = withing group variance; tau00 = between group variance; ICC = proportion of variance explained for each group level


#Sum of squares to give the percentage explained by each variable in the model####
#Day
global_mdl_inter_day_Seiche_aov <-  anova(global_mdl_inter_day_seiche)

#Percentage explained by each variable
round(global_mdl_inter_day_Seiche_aov$`Sum Sq`[1]*100/sum(global_mdl_inter_day_Seiche_aov$`Sum Sq`),2)
round(global_mdl_inter_day_Seiche_aov$`Sum Sq`[2]*100/sum(global_mdl_inter_day_Seiche_aov$`Sum Sq`),2)
round(global_mdl_inter_day_Seiche_aov$`Sum Sq`[3]*100/sum(global_mdl_inter_day_Seiche_aov$`Sum Sq`),2)
round(global_mdl_inter_day_Seiche_aov$`Sum Sq`[4]*100/sum(global_mdl_inter_day_Seiche_aov$`Sum Sq`),2)

#Night
global_mdl_inter_night_seiche_aov <-  anova(global_mdl_inter_night_seiche)

#Percentage explained by each variable
round(global_mdl_inter_night_seiche_aov$`Sum Sq`[1]*100/sum(global_mdl_inter_night_seiche_aov$`Sum Sq`),2)
round(global_mdl_inter_night_seiche_aov$`Sum Sq`[2]*100/sum(global_mdl_inter_night_seiche_aov$`Sum Sq`),2)
round(global_mdl_inter_night_seiche_aov$`Sum Sq`[3]*100/sum(global_mdl_inter_night_seiche_aov$`Sum Sq`),2)
round(global_mdl_inter_night_seiche_aov$`Sum Sq`[4]*100/sum(global_mdl_inter_night_seiche_aov$`Sum Sq`),2)





#Random intercept per fishid with interactions
tic("model fitting")
global_mdl2 <- lmer(formula = det_depth ~ diel_period * (lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center) + (1|fishid),
                    data = detections_mdl, REML = T, 
                    lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #random group intercept
toc()



AICctab(global_mdl1, global_mdl1_thick, global_mdl2)

#Random slope per fishid
tic("model fitting")
global_mdl3 <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center + (1 + lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center | fishid),
                    data = detections_mdl, REML = T,
                    lmerControl(optimizer = 'bobyqa')) #random group slope
toc()



#random intercept and random slope per fish id - Day (Seiche only) #####
tic('Fuck, it took this amount of time to run the model')
mdl_random_slo_int_day <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1+lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center|fishid),
                               data = detections_mdl[diel_period=='day' & Seiche == 'on'], REML = T, 
                               lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()

#random intercept and random slope per fish id - Night (Seiche only) #####
tic('Fuck, it took this amount of time to run the model')
mdl_random_slo_int_night <- lmer(formula = det_depth ~ lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center + (1+lake_therm_thickness_smoothed + therm_temp_strength + lake_therm_depth_smoothed_center + det_therm_deviation_center|fishid),
                                 data = detections_mdl[diel_period=='night' & Seiche == 'on'], REML = T, 
                                 lmerControl(optimizer = 'bobyqa', calc.derivs = FALSE)) #uncorrelated random intercept and random slope within fishid
toc()


tab_model(mdl_random_slo_int_day, mdl_random_slo_int_night, dv.labels = c("Day","Night"), 
          title='Wels (Seiche only)', p.style='asterisk',
          pred.labels = c('Intercept', 'Thermocline thickness (m)',
                          'Thermocline strength (ºC)',
                          'Thermocline depth (m)',
                          'Seiche strength'))#sigma² = withing group variance; tau00 = between group variance; ICC = proportion of variance explained for each group level





#Sum of squares to give the percentage explained by each variable in the model####
#Day
mdl_random_slo_int_day_aov <-  anova(mdl_random_slo_int_day)



#Percentage explained by each variable
round(mdl_random_slo_int_day_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_day_aov$`Sum Sq`),2)
round(mdl_random_slo_int_day_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_day_aov$`Sum Sq`),2)
round(mdl_random_slo_int_day_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_day_aov$`Sum Sq`),2)
round(mdl_random_slo_int_day_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_day_aov$`Sum Sq`),2)

#Night
mdl_random_slo_int_night_aov <-  anova(mdl_random_slo_int_night)

#Percentage explained by each variable
round(mdl_random_slo_int_night_aov$`Sum Sq`[1]*100/sum(mdl_random_slo_int_night_aov$`Sum Sq`),2)
round(mdl_random_slo_int_night_aov$`Sum Sq`[2]*100/sum(mdl_random_slo_int_night_aov$`Sum Sq`),2)
round(mdl_random_slo_int_night_aov$`Sum Sq`[3]*100/sum(mdl_random_slo_int_night_aov$`Sum Sq`),2)
round(mdl_random_slo_int_night_aov$`Sum Sq`[4]*100/sum(mdl_random_slo_int_night_aov$`Sum Sq`),2)

plot_model(mdl_random_slo_int_day, title = 'Wels (day)')
plot_model(mdl_random_slo_int_night, title = 'Wels (night)')

plot_model(mdl_random_slo_int_night, type = "pred", axis.lim = c(1,5),
           terms = c('lake_therm_depth_smoothed_center', 'fishid [T412073, T412074, T412075, T412076]'))

plot_residuals(mdl_random_slo_int_night)




#Running more than 1 model at once at different cores
models <- c('lmer(formula = det_depth ~ lake_therm_thickness_smoothed + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center + (1|fishid), data = i, REML = F)',
            'lmer(formula = det_depth ~ lake_therm_thickness_smoothed + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center + (1 + lake_therm_thickness_smoothed + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center | fishid), data = i, REML = F)')
m.list <- f_lmer_mc(data = detections_mdl, calls = models, mc.cores = 20)


#Diagnostic plot
sjp.lmer(global_mdl8, type = 'fe')
sjp.lmer(global_mdl_inter_night_seiche_aov, type = "fe.ri", sort='fishid')

plot_model(global_mdl_inter_night_seiche, type = 'est', title = 'Tench', vline.color = "black")

#Model assumptions checking####
#homoscedasticity test
car::leveneTest(residuals(global_mdl6) ~ detections_mdl$lake_therm_thickness_smoothed + detections_mdl$lake_therm_depth_smoothed_center + detections_mdl$diel_period + detections_mdl$det_therm_deviation_center)

#Normality of residuals
qqnorm(residuals(global_mdl6))#OK

#Linearity in each variable
#Check for each one of the independent variables
##Thermocline thickness
ggplot(data.frame(x1= detections_mdl$lake_therm_thickness_smoothed,
                  pearson = residuals(global_mdl6,type="pearson")), aes(x = x1, y = pearson))+
  geom_point(alpha=0.5, shape='.')+
  geom_smooth(method='lm', se=F)+
  theme_bw()
##Balanced thermocline depth (center)
ggplot(data.frame(x1= detections_mdl$lake_therm_depth_smoothed_center,
                  pearson = residuals(global_mdl6,type="pearson")), aes(x = x1, y = pearson))+
  geom_point(alpha=0.5, shape='.')+
  geom_smooth(method='lm', se=F)+
  theme_bw()
##Period of the day
ggplot(data.frame(x1= detections_mdl$diel_period,
                  pearson = residuals(global_mdl6,type="pearson")), aes(x = x1, y = pearson))+
  geom_point(alpha=0.5, shape='.')+
  geom_smooth(method='lm', se=F)+
  theme_bw()
##Seiche strength
ggplot(data.frame(x1= detections_mdl$det_therm_deviation_center,
                  pearson = residuals(global_mdl6,type="pearson")), aes(x = x1, y = pearson))+
  geom_point(alpha=0.5, shape='.')+
  geom_smooth(method='lm', se=F)+
  theme_bw()


#Selection of the best model
AICtab(global_mdl1, global_mdl2, global_mdl3, global_mdl4)#Akaike Information Criterion
BICtab(global_mdl1, global_mdl2, global_mdl3, global_mdl4)#Bayesian Information Criterion


#Inter Class Correlation (ICC) - (Measure of how much of the variation in the response variable is accounted by the random effect)
r1Var <- as.numeric(VarCorr(global_mdl1)[["fishid"]])
residVar <- attr(VarCorr(global_mdl1), "sc")^2
round((r1Var/(r1Var+residVar)), 2) #percentage of variation explained by the random effect

AICtab(global_mdl8, global_mdl6)

#Model diagnosis
ranef(global_mdl1)#estimated deviation between each level of the random effect and the overall average:
fixef(global_mdl1)#fixed effects outputs
round(MuMIn::r.squaredGLMM(global_mdl1),2)#R2m = R2 for the marginal test (fixed); R2c = R2 for the conditional test (fixed + random)

# #Model simplification by backwards
# lmerTest::step(model=global_mdl2)#the analysis suggested to keep all variables into the explanatory side of the model 
# #test the significance of each variable (full model x reduced model)
# #Temperature
# asp_mdl_int_temp <- lmer(formula = eggs_sca ~ dayoy_sca + visit_sca + dayNight + (1|year),
#                          data = df_complete, REML = F) #intercept varying among years and among fish_id within years (nested random effects)
# anova(asp_mdl_int,asp_mdl_int_temp)#Temperature is significant
# #Period of the day
# asp_mdl_int_day <- lmer(formula = eggs_sca ~ dayoy_sca + visit_sca + temp_sca + (1|year),
#                         data = df_complete, REML = F) #intercept varying among years and among fish_id within years (nested random effects)
# anova(asp_mdl_int,asp_mdl_int_day)#Period of the day is significant
# #Bleak presence
# asp_mdl_int_bleak <- lmer(formula = eggs_sca ~ dayoy_sca + dayNight + temp_sca + (1|year),
#                           data = df_complete, REML = F) #intercept varying among years and among fish_id within years (nested random effects)
# anova(asp_mdl_int,asp_mdl_int_bleak)#Bleak presence is significant
# #day of the year
# asp_mdl_int_dayoy <- lmer(formula = eggs_sca ~ visit_sca + dayNight + temp_sca + (1|year),
#                           data = df_complete, REML = F) #intercept varying among years and among fish_id within years (nested random effects)
# 
# anova(asp_mdl_int, asp_mdl_int_dayoy)#Day of the year is significant

#stepwise selection
MASS::stepAIC(object = global_mdl6, direction = 'backward')

#Models with each individual in separate####
tench_ids <- unique(detections_mdl[,fishid])#vector of unique fish ids

mdl_list <- list()#creating a list to store the models

#for cycle for the computation of all models 
for (i  in 1:length(unique(detections_mdl[,fishid]))) {
  mdl_list[[i]] <- lm(formula = det_depth ~ therm_temp_strength + lake_therm_thickness_smoothed + lake_therm_depth_smoothed_center + diel_period + det_therm_deviation_center, data = detections_mdl[fishid==tench_ids[[i]]])
}

#apllying the summary of all models at once
lapply(mdl_list , summary)

#creating an automatic export of the results into a table
sjt.lm(mdl_list[[1]],mdl_list[[2]],mdl_list[[3]],mdl_list[[4]],mdl_list[[5]],
       mdl_list[[6]],mdl_list[[7]],mdl_list[[8]],mdl_list[[9]],mdl_list[[10]],
       mdl_list[[11]],mdl_list[[12]],mdl_list[[13]])