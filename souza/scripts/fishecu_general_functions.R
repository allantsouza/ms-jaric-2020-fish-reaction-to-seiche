#Handy funcitons used for all converting of samplings for FISHECU


#fix name of column, print error if there are more than 1 column with pattern specified
fixColumnName <- function(x, xpattern, new.name){
  colidx <- grep(pattern = xpattern, x = names(x)); 
  if(! new.name %in% names(x)){
    if(length(colidx) > 1) print(paste("ERROR:more columns starting with pattern -", xpattern));
    if(length(colidx) == 0){x[,new.name] <- NA}else{names(x)[colidx] <- new.name};
  }
  return(x)  
}

#fix species
fixSpecies <- function(x, species.typos, species.db){
  x$Species <- sub("\\\u009d$", "", x$Species)
  #load map of species - first correct typos
  x <- merge(x, species.typos, by  = "Species", all.x = T)
  x[!is.na(correctsp), Species := correctsp]
  sort(unique(x$Species))
  species.unq <- sort(unique(x$Species))
  x <- merge(x, species.db[, .(czechcodename, czechname)], by.x = c("Species"), by.y= c("czechcodename"), all.x = T)
  #CHECKs
  #get species which are not in species table
  print(paste("Following species were not found in database:" , paste(unique(x[is.na(czechname)]$Species), collapse = ","))) # '0' OK
  x
}

#function to get for each input time if it was night or day
#input:
# x - vector of type POSIXct
# lat - numeric vector of lenght one
# lon - numeric vector of lenght one
# sunsetsunrise <- suncalc::getSunlightTimes(date = as.Date(c("2017-05-02", "2017-05-05", "2017-05-06")), lat = 35, lon = 32, keep = c("sunrise", "sunset"))[,c(4,5)]
# sunrise.sunset.t <- data.table(sunrise.set(49+(34/60)+((42/60)/60), 15+(15/60)+((21/60)/60), "2017-01-01" , num.days = 365, timezone = "UTC"))
# sunrise.sunset.t[, rize.sec := as.numeric(difftime(sunrise, round.POSIXt(x = sunrise, "days"), units = "secs"))]
# sunrise.sunset.t[, set.sec := as.numeric(difftime(round.POSIXt(x = sunset, "days"), sunset, units = "secs"))]
# 
# thresholds <- sunrise.sunset.t[, .(first.date = min(as.Date(sunrise)), minrise = min(rize.sec), maxrise = max(rize.sec), minset = min(set.sec), maxset = max(set.sec)), by = .(month = month(sunrise))] 
# thresholds[, data.table(period = rep(c("eveningshift", "night", "morningshift"), each = 4), 
#                         x = c(minset, minset, maxset, maxset,
#                               maxset, maxset, minrise, minrise, 
#                               minrise, minrise, maxrise, maxrise),
#                         y = rep(c(-Inf, Inf, Inf, -Inf), 3)), by = .(first.date)]
# 



getNightDay <- function(x, lat = 49.5765639, lon = 14.6637706, xtz = "UTC", returnSunrSuns = F, day = "day", night = "night"){
  if(length(x) == 0) stop("length of input is 0")
  if(all(is.na(x))) return(as.character(x))
  
  library(StreamMetabolism)
  library(data.table)
  from <- min(x, na.rm = T)
  to <- max(x, na.rm = T)
  sunrise.sunset <- as.data.table(sunrise.set(lat, lon, from , num.days = difftime(time1 = to, time2 = from, units = "days", tz = xtz)+2, timezone = xtz))
  sunrise.sunset$Date <- as.Date(sunrise.sunset$sunrise, tz = xtz)
  res <- merge(data.table(x, xorder = 1:length(x), Date = as.Date(x)), sunrise.sunset, by = "Date", all.x = T)
  res[, Diel.period := as.character(ifelse(x > sunrise & x < sunset, "day", "night"))]
  setkey(res, xorder)
  if(returnSunrSuns){
    return(res[, .(Diel.period, sunrise, sunset)])
  }else{
    return(res$Diel.period)  
  }
}



#recursive function to get samplings close to each other
getOvelappingSamplings<- function(samplingids, dates, days.distance = 10){
  library(data.table)
  if(length(samplingids) == 0 | length(samplingids) == 0) {print("ERROR: vector of length 0 as input"); return(NA)}
  x <- data.table(samplingid = samplingids, Date = dates)
  x$processed <- F
  x$group <- -1
  group.i <- 1
  for(i in 1:nrow(x)){
    if(x$processed[i] == F){
      oversamplings.idxs <- getOverSamplingsRec(i, x, days.distance)  
      x[oversamplings.idxs, processed := T]
      x[oversamplings.idxs, group := group.i]
      group.i <- group.i + 1
    }  
  }
  return(x[, .(group)])
}

getOverSamplingsRec <- function(idx, aa, days.distance){
  aa$processed[idx] <- T
  idxs <- unique(which(aa$processed == F & aa$Date %between% c(aa[idx]$Date-days.distance, aa[idx]$Date+days.distance)))
  if(length(idxs) == 0) return(idx)
  aa$processed[idxs] <- T
  return(unique(c(idx, unlist( sapply(unique(idxs), function(idx, aa, days.distance) getOverSamplingsRec(idx, aa, days.distance), days.distance = days.distance, aa = aa, simplify = T), recursive = TRUE))))
}

#Fix chzech letters - substitute special czech letters
#@input: character vector to be cleaned
#@output: cleaned character vector
fixCzechLetters <- function(x){
  #Vowels
  x <- gsub(pattern = "\u00E1", replacement = "a", x )#Lowercase long a (á)
  x <- gsub(pattern = "\u00C1", replacement = "A", x )#Uppercase long a (Á)
  x <- gsub(pattern = "\u00E9", replacement = "e", x )#Lowercase long e (é)
  x <- gsub(pattern = "\u00C9", replacement = "E", x )#Uppercase long e (É)
  x <- gsub(pattern = "\u011B", replacement = "e", x )#Lowercase soft e (ě)
  x <- gsub(pattern = "\u011A", replacement = "E", x )#Uppercase soft e (Ě)
  x <- gsub(pattern = "\u00ED", replacement = "i", x )#Lowercase long i (í)
  x <- gsub(pattern = "\u00CD", replacement = "I", x )#Uppercase long i (Í)
  x <- gsub(pattern = "\u00F3", replacement = "o", x )#Lowercase long o (ó)
  x <- gsub(pattern = "\u00D3", replacement = "O", x )#Uppercase long o (Ó)
  x <- gsub(pattern = "\u00FA", replacement = "u", x )#Lowercase long u (ú)
  x <- gsub(pattern = "\u00DA", replacement = "U", x )#Uppercase long u (Ú)
  x <- gsub(pattern = "\u016F", replacement = "u", x )#Lowercase ring u (ů)
  x <- gsub(pattern = "\u016E", replacement = "U", x )#Uppercase ring u (Ů)
  #Consonants
  x <- gsub(pattern = "\u010D", replacement = "c", x )#Lowercase soft c (č)
  x <- gsub(pattern = "\u010C", replacement = "C", x )#Uppercase soft c (Č)
  x <- gsub(pattern = "\u010F", replacement = "d", x )#Lowercase soft d (ď)
  x <- gsub(pattern = "\u010E", replacement = "D", x )#Uppercase soft d (Ď)
  x <- gsub(pattern = "\u0148", replacement = "n", x )#Lowercase soft n (ň)
  x <- gsub(pattern = "\u0147", replacement = "N", x )#Uppercase soft n (Ň)
  x <- gsub(pattern = "\u0159", replacement = "r", x )#Lowercase soft r (ř)
  x <- gsub(pattern = "\u0158", replacement = "R", x )#Uppercase soft r (Ř)
  x <- gsub(pattern = "\u0161", replacement = "s", x )#Lowercase soft s (š)
  x <- gsub(pattern = "\u0160", replacement = "S", x )#Uppercase soft s (Š)
  x <- gsub(pattern = "\u0165", replacement = "t", x )#Lowercase soft t (ť)
  x <- gsub(pattern = "\u0164", replacement = "T", x )#Uppercase soft t (Ť)
  x <- gsub(pattern = "\u00FD", replacement = "y", x )#Lowercase long y (ý)
  x <- gsub(pattern = "\u00DD", replacement = "Y", x )#Uppercase long y (Ý)
  x <- gsub(pattern = "\u017E", replacement = "z", x )#Lowercase soft z (ž)
  x <- gsub(pattern = "\u017D", replacement = "Z", x )#Uppercase soft z (Ž)
  return(x)
}

#Function to detect outliers in LW relationship  
detectOuts <- function(W, L, rlm.thresh){
  library(MASS)
  if(length(W) >1 & length(L) > 1){
    m.rlm <- rlm(log10(W) ~ log10(L), maxit = 200)
    return(m.rlm$w < rlm.thresh)
  }else{return(F)}
}

#general function computing value per unit of effort
#input - dt of samplings, catch, vector of catch split factors, vector of sampling split factors, samplingid column, value.var is name of column containing count/weight/whatever

getVPUE <- function(samplings, catch, split.factors.catch, split.factors.samplings, id.colname = "samplingid", effort.colname = "Effort", returnlist = F, value.var){
  #check input
  library(data.table)
  if(!is.character(split.factors.catch) | !is.character(split.factors.samplings)){stop("Split.factors argument is expected to be character vector")  }
  if(!(is.data.table(samplings) & is.data.table(catch))){stop("Catch and sampling has to be data.table object")}
  if(length(id.colname) != 1){stop("Wrong length of id.colname argument: expected length == 1")}
  if(!value.var %in% names(catch)){stop(paste("The value.var:", value.var, " column is missing in catch argument"))}
  if(!id.colname %in% names(samplings)){stop(paste("The id.colname:", id.colname, " column is missing in samplings argument"))}
  if(!id.colname %in% names(catch)){stop(paste("The id.colname:", id.colname, " column is missing in catch argument"))}  
  if(!effort.colname %in%  names(samplings)){stop(paste("The effort.colname:", effort.colname, " column is missing in samplings argument"))}
  
  #check if there are some fish
  if(nrow(catch) == 0){warning("There was no fish caught (this might mean that sampling ids are not matching!)")}
  
  #dataframes should 
  
  #TODO: check if all split factors are present in tables
  #add samplingid column for easy reference
  catch$samplingid <- catch[, id.colname, with = F]
  samplings$samplingid <- samplings[, id.colname, with = F]
  samplings$Effort <- samplings[, effort.colname, with = F]
  catch$value.var <- catch[, value.var, with = F]
  
  #subset catch only to catch from sampling table
  catch <- catch[samplingid %in% samplings$samplingid]
  
  #split into one vector
  split.factors <- c(split.factors.catch, split.factors.samplings)
  
  #keep only columns needed for computation
  catch <- catch[, c("samplingid", "value.var", split.factors.catch), with = F]
  samplings <- samplings[, c("samplingid","Effort", split.factors.samplings), with = F]
  
  #add number of samplings in factor group
  samplings[, sampcount := .N, by = split.factors.samplings]
  
  #if there are no catch factors 
  if(length(split.factors.catch) > 0){
    catch.uniques <- do.call(what = CJ, args =  c(catch[, split.factors.catch, with = F], list(unique = T)))
    #add all species
    samplings.sp <- CJ.tables(samplings, catch.uniques)
  }else{ samplings.sp <- samplings}
  
  #if all fish dont have weight - return error
  catch.sum <- catch[, .(value.var.sum = sum(value.var), fishcount = .N), by = c("samplingid", split.factors.catch)]
  catch.sampl <- merge(catch.sum, samplings.sp, all = T, by = names(samplings.sp)[names(samplings.sp) %in% names(catch.sum)])
  catch.sampl[is.na(value.var.sum), value.var.sum :=  0]
  catch.sampl[is.na(fishcount), fishcount :=  0]
  catch.sampl[, VPUE := value.var.sum / Effort]
  
  #populate sampling table by factor values
  catch.agg <- catch.sampl[, .(VPUE.mean = mean(VPUE), VPUE.sd = sd(VPUE), VPUE.sum = sum(value.var.sum),  sampcount = sampcount[1], catch.rowcount = sum(fishcount)), by = split.factors]
  
  names(catch.agg) <- gsub(x = names(catch.agg), pattern = "VPUE", replacement = value.var)
  #0 where the sampling was done
  #NA where the sampling was not done
  if(returnlist == T){
    return(list(VPUE.per.net = catch.sampl[, c("samplingid", split.factors, "VPUE"), with = F], VPUE.per.factors  = catch.agg ))
  }else{
    return(catch.agg)
  }
}



#function to cross join two data.tables
CJ.tables <- function(X,Y){ 
  unique_name <- last(make.unique(c(colnames(X),colnames(Y),"temporalColumnForCJtables"))) 
  return(X[,c(setNames(1,unique_name),.SD)][Y[,c(setNames(1,unique_name),.SD)],on=unique_name,allow.cartesian=TRUE]) 
}
#round posixct time to specific bin in seconds
round_POSIXct <- function(x, bin.secs){
  if(class(x)[1] != "POSIXct") stop("parameter x has to be type POSIXct")
  library(lubridate)
  as.POSIXct(round(as.numeric(x)/bin.secs)*bin.secs,origin= "1970-01-01", tz = tz(x))
}
floor_POSIXct <- function(x, bin.secs){
  if(class(x)[1] != "POSIXct") stop("parameter x has to be type POSIXct")
  library(lubridate)
  as.POSIXct(floor(as.numeric(x)/bin.secs)*bin.secs,origin= "1970-01-01", tz = tz(x))
}
ceiling_POSIXct <- function(x, bin.secs){
  if(class(x)[1] != "POSIXct") stop("parameter x has to be type POSIXct")
  library(lubridate)
  as.POSIXct(ceiling(as.numeric(x)/bin.secs)*bin.secs,origin= "1970-01-01", tz = tz(x))
}

#Functionto merge 12 and 4 size mesh panel into 16
#load 12
#load 4
#check if stratum, locality is the same for the pairs, 
#check if 
#divide 4 by 4
#union fish
#take date from 12
#leave GPS, notes, depthmin, depthmax etc. as it is. Its up to every one to 


#Shanon index
#Function computing shanon index 
#input: vector defining a species, vector defining abundances of each species(optional)
#if the abundance is not supplied, abundance is computed as number of occurences of each species in vector species
#NOTE: one species can occure in species vector more than once (can be used with pure catch table)
#output:Shanon index
#Example:
#catch[, .(Shannon.idx = computeShannon(species, abundance), by = .(reservoir, locality)]

computeShannon <- function(species,  abundance = NULL){
  if(is.null(abundance)){
    x <- data.table(species = species)
    x.sp <- x[, .(spcount = .N), by = species]
  }else{
    x <- data.table(species = species, abundance = abundance)
    x.sp <- x[, .(spcount = sum(abundance)), by = species]
  }
  x.sp[, totcount := sum(spcount)]
  x.sp[,  p.sp := spcount / totcount]
  return(x.sp[, -sum(p.sp*log(p.sp))])
}

#roll join function by constant width time window
#input: x - vector of values to be smoothed
#time.vec - vector of timestamps
#win.widt - width of the window in seconds
#FUN - function to apply on values in the window
rollTimeWindow <- function(x, time.vec, win.width, FUN){
  x.out <- vector(mode = class(x), length = length(x))
  x.out[1:length(x.out)] <- NA
  if(is.unsorted(time.vec)){
    stop("Time vector is not in increasing order!")
  }
  FUN <- match.fun(FUN)
  time.vec.minus <- time.vec-win.width
  time.vec.plus <- time.vec+win.width
  for(i in 1:length(x.out)){
    x.out[i] <- FUN(x[which(time.vec > time.vec.minus[i] & time.vec < time.vec.plus[i])])
  }
  return(x.out)
}


get1dKernel <- function(x, lev, bw, from, to, getPlot = F, return_npatches = F){
  if(any(x < from | x > to)) stop("Some samples do not fall into specified boundaries range (from, to)")
  #reflect dataset both sides behind boundaries - mirror only points in distance equal to 10% of the interval size
  mirror.dist <- (to-from) * .05
  x2 <- c(x,  from - x[x - from  < mirror.dist ], to + (to - x[to - x  < mirror.dist]))
  dens.f <- density(x2, kernel= "gaussian", bw = bw, from = from, to = to)
  if(getPlot){
    hist(x2, xlim= c(from, to), col = "red", breaks = (to-from)/250)
    par( new = T)
    plot(dens.f, lwd = 3)
  }
  density_bars_sorted <- sort(dens.f$y, decreasing = FALSE)
  y.cumsum <- cumsum(density_bars_sorted)
  #since the mirroring is used (almost twice as much points on negative side), the percentage out has to be devided by 2
  accept.cum.level <- max(y.cumsum)*((100-lev)/100)
  #get number of bins of density function where the probability is more than percentage level
  kernel.bins <- length(y.cumsum[y.cumsum > accept.cum.level])
  #multiply by real size of bin in meters
  kernel_width <- kernel.bins*(dens.f$x[2]-dens.f$x[1])
  
  if(return_npatches){
    accept.level <- density_bars_sorted[length(y.cumsum)-kernel.bins]
    isPeak <- dens.f$y > accept.level
    isPeak_idx <- rleid(isPeak)
    #in case when the first group is peak
    if(isPeak[1]) isPeak_idx <- isPeak_idx + 1
    n_patches <- max(isPeak_idx[isPeak_idx %%2 == 0]/2)
    plot(dens.f$x, dens.f$y, col= (dens.f$y > accept.level)+1, pch = 20)
    abline(h= accept.level)
    return(kernel_width, n_patches)
  }else{
    return(kernel_width)
  }
}



getNighttimePols <- function(sunsetsunrise){
  nightime.vector <- sort(c(sunsetsunrise[-1,1], sunsetsunrise[-nrow(sunsetsunrise),2]))
  x <- rep(nightime.vector, each = 2)
  y <- rep(c(-Inf,Inf), times = length(nightime.vector))
  xord <- rep(c(1,2,4,3), times =length(nightime.vector)/2)
  night.polygons <- data.table(x,y, id = rep(1:(round(length(nightime.vector)/2)), each = 4), xord = xord)
  setkey(night.polygons, id,xord)
  return(night.polygons)
}