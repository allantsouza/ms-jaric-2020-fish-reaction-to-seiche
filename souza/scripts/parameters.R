#Parameters
step_order_par <- 1
degree_per_meter_par <- 2
DateRange <- range("2015-06-16","2015-10-15")

#Specify species and lake
  #specify 1 species in question (tench, pike, wels, rudd, roach, perch)
  species <- c('wels')#tench
  #specify 1 lake in question (Chabarovice, Most)
  lake <- "Chabarovice"
  
#Far merge exclusion
  #parameter used for excluding detections which are further than X seconds from the closest position
  excl_dets_further_than_pos <- 60*15
  #exclude detections for which is the thermocline detected futher than X seconds
  excl_dets_further_than_therm <- 60*15

#Balance depth of thermocline
  #rolling window width for smoothing thermocline balance
  therm_bal_ws <- 60*60*24*4