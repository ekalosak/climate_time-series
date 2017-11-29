# Author: Eric Kalosa-Kenyon
# Load data

# rm(list=ls()) # clean workspace

library(ncdf4) # load .nc data files

# Setup parameters
data_dir = "~/Desktop/STA237A/Project/"
#data_dir = "~/Documents/Classes/STA237A/project/data/"
data_file = "lai_Lmon_CCSM4_rcp45_r2i1p1_200601-210012.nc"
data_path = paste(data_dir, data_file, sep="")
attr_name = "lai" # leaf area index
plot.bool = T # plot results along the way?
using.rstudio = F

#Setwd to current source location
if(using.rstudio){
    src_dir=dirname(rstudioapi::getActiveDocumentContext()$path)
}else{
    src_dir="./"
}
setwd(src_dir)

# Load the.nc
ncin = nc_open(data_path)

# Extract lat, lon, time
lon = ncvar_get(ncin, "lon")
lat = ncvar_get(ncin, "lat")
loc.af = c(
           min(which(lat > 3)),
           min(which(lon > 25))
           )
time.mo = ncvar_get(ncin, "time")
tunits = ncatt_get(ncin, "time", "units")

# Extract attribute and data
attr = ncatt_get(ncin, attr_name)
attr_data = ncvar_get(ncin, attr_name)
d = attr_data # shorthand for the array of space X time
# NOTE: d is indexed d[lon, lat, time], but the natural location indexing is
# loc = [lat, lon] so to index d by a location you must do
# d[loc[2], loc[1], time] rather than d[loc[1], loc[2], time]

## Extract the data for the time series in a single location in sub-saharan
## Africa (i.e. single v=LAI, single s=African, all T=2000 to 2100)
d.af = d[loc.af[2],loc.af[1],] # d.af is african single loc data from 2000 to 2100

# SOURCES:
# 1. http://geog.uoregon.edu/GeogR/topics/netCDF-read-ncdf4.html
