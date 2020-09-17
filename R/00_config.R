# Load libraries
library(renv)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(lubridate)
library(RPostgreSQL)
library(data.table)
library(here)


# Global variables -------------------------------------------------------------

DATE_RANGE <- c("2015-07-05","2015-07-25")

PAR_THERMOCLINE_SLOPE <- 1

#Identification of fast Seiche shift periods - "shift speed" has to be larger than a threshold speed (PAR_SHIFT_HEIGHT, which is X meters per hour)
PAR_SHIFT_SPEED <- 0.5

#Balance depth of thermocline
#rolling window width for smoothing thermocline balance
PAR_THERMOCLINE_BALACE <- 60*60*24*4

# TODO connection to database
con <- dbConnect(drv = PostgreSQL(), 
                 dbname = Sys.getenv("dbname"),
                 host= Sys.getenv("dbhost"), 
                 user= Sys.getenv("dbuser"),
                 password = Sys.getenv("dbpassword"))


# global ggplot theme
theme_set(theme_minimal())


# Detections --------------------------------------------------------------

# All positions are projected on lake axis but to keep as much detections as possible, take also detections which had closest position X seconds away
PAR_DET_POS_TIMEDIFF <- 30*60
