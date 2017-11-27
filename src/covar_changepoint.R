# Author: Eric
# Calculate location-covariance matrices over time
# and extract a time-series signal from them

rm(list=ls())
source("detrend.R") #, echo=TRUE)

lat = dim(d)[1] # TODO: check that lat=1 and lon=2, or reverse?
lon = dim(d)[2]
obs = dim(d)[3] # number of observations in data frame
yr = 12 # observations are every month, so 12 of them per year
yrs = obs/yr # number of years simulated

covts = array(NA, dim=c(lat*lon, lat*lon, yrs))

# TODO: CONTINUE: For each year of observations, calculate locationwise
# covariance
