# Author: Eric
# Calculate location-covariance matrices over time
# and extract a time-series signal from them

rm(list=ls())
source("detrend.R") #, echo=TRUE)

nlat = dim(d)[1] # TODO: check that lat=1 and lon=2, or reverse?
nlon = dim(d)[2]
nobs = dim(d)[3] # number of observations in data frame
yr = 12 # observations are every month, so 12 of them per year
nyrs = nobs/yr # number of years simulated

# Construct the covariance time series array.
# each slice in the 2rd dimension is the location covariance of LAI in the year
# sliced. covts uses column-major indexing i.e. cv_(i,j,k) = cv[i + n*j, k]
# covts = Matrix(0, dim=c((lat*lon)**2, yrs), sparse=TRUE)
# for(y in 1:yrs){
# }

# TODO: CONTINUE: For each year of observations, calculate locationwise
# covariance

# TODO: create ST data for the spatiotemporal package
library(SpatioTemporal)
loc_ids = 1:(nlat*nlon) # names of each grid cell
covars = data.frame(ID=loc_ids)
# rownames(covars) = covars$ID
covars[,c("lat", "long")] = expand.grid(lat, lon) # all pairs of lat/lon
obs = data.frame(row.names=time.mo) # ID by time containing observations
