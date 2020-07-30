#Libraries
library(data.table)
library(ggplot2)
library(ggthemes) # extension of ggplots, used for home range visualsation
library(lubridate)
library(RPostgreSQL)
library(StreamMetabolism) # for calculation of sunset and sunrise 
library(shiny)
source(file = 'scripts/fishecu_general_functions.R', local = T)
library(nlme)
library(MASS)
#install.packages('shiny')
