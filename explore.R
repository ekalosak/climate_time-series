# Author: Eric Kalosa-Kenyon
#
# Opens and plots some leaf area index data

rm(list=ls()) # clean workspace

library(ncdf4) # load .nc data files

# Setup parameters
data_dir = "./data/"
data_file = "lai_Lmon_CCSM4_rcp45_r2i1p1_200601-210012.nc"
data_path = paste(data_dir, data_file, sep="")
attr_name = "lai" # leaf area index
plot.bool = T # plot results along the way?
plot.filename = "basic_analysis.pdf"

# Open outfile in write mode
pdf(file=plot.filename)

# Load the.nc
ncin = nc_open(data_path)

# Extract lat, lon, time
lon = ncvar_get(ncin, "lon")
lat = ncvar_get(ncin, "lat")
loc.af = c(
           min(which(lat > 3)),
           min(which(lon > 25))
           )
t = ncvar_get(ncin, "time")
tunits = ncatt_get(ncin, "time", "units")

# Extract attribute and data
attr = ncatt_get(ncin, attr_name)
attr_data = ncvar_get(ncin, attr_name)
d = attr_data # shorthand for the array of space X time
# NOTE: d is indexed d[lon, lat, time], but the natural location indexing is
# loc = [lat, lon] so to index d by a location you must do
# d[loc[2], loc[1], time] rather than d[loc[1], loc[2], time]

## Plot the global LAI at time 0
# extract variables to plot
x = lon # 288 longitudes
y = lat # 192 latitudes
z = d[,,1] # all LAI values across the planet at time 1

# plotting subroutines
leaf.colors = function(x){rev(terrain.colors(x))}
if(plot.bool){
    filled.contour(x,y,z, color = leaf.colors, asp = 1,
                   plot.axes={
                       axis(1); axis(2);
                       points(
                              floor(lon[loc.af[2]]),
                              floor(lat[loc.af[1]])
                              )
                   })
    title(xlab="Longitude", ylab="Latitude", main="LAI at T=1")
}

## Extract the data for the time series in a single location in sub-saharan
## Africa (i.e. single v=LAI, single s=African, all T=2000 to 2100)
d.af = d[loc.af[2],loc.af[1],] # d.af is african single loc data from 2000 to 2100
if(plot.bool){
    plot(t, d.af, xlab="Time", ylab="Tree cover (LAI)",
          main="Single location sub-Saharan African tree cover over a century",
          type="l")
}

## Extract a single day each year as an attempt to de-season the data
st.day = 1 # starting day, 1 is Jan 01, 2005
# NOTE: the following doesn't work... need to go throught the "t" vector that
# holds the times of each observation
day.ixs = seq(from=st.day, to=length(d.af), by=365)

## Standardize variance, detrend, deseason etc. standard time series analysis

# Check acf and pacf
if(plot.bool){
    par(mfrow=c(2,1))
    acf(d.af, main="Raw data (ACF)")
    pacf(d.af, main="Raw data (PACF)")
}

# detrend with polynomial model
df = data.frame(obs=d.af, time=t)
degree = 6
fit.lm = lm(obs ~ poly(time, degree, raw=TRUE), df)
d.af.detr = d.af-predict(fit.lm)

if(plot.bool){
    par(mfrow=c(3,1))
    plot(t, predict(fit.lm),
         xlab="Time", ylab="Polynomial signal",
         type="l")
    title(main="Polynomial trend")

    plot(t, d.af,
         xlab="Time", ylab="True signal",
         type="l", col="red")
    title(main="Observations")

    plot(t, d.af.detr,
         xlab="Time", ylab="Detrended signal",
         type="l", col="blue")
    title(main="Detrended")
    par(mfrow=c(1,1))
}

# Check acf and pacf with detrended data
if(plot.bool){
    par(mfrow=c(2,1))
    acf(d.af.detr, main="Detrended data (ACF)")
    pacf(d.af.detr, main="Detrended data (PACF)")
}

# try using differences
d.af.diff = diff(d.af)
if(plot.bool){
    par(mfrow=c(2,1))
    acf(d.af.diff, main="Differenced data (ACF)")
    pacf(d.af.diff, main="Differenced data (PACF)")
}

# moving average to de-season
filter.length = 12 # 12 used because t[1]-t[13] = -365
filt.ma = filter(d.af.detr, rep(1,filter.length)/filter.length)
if(plot.bool){
    par(mfrow=c(1,1))
    plot(t, d.af.detr, type="l", xlab="Days since 01-01-2005", ylab="LAI")
    lines(t, filt.ma, type="l", col="red")
    title(main="One year seasonality removal from detrended observations")
}

d.af.deseas = d.af.detr - filt.ma
if(plot.bool){
    par(mfrow=c(3,1))
    plot(t, d.af.deseas,
         type="l", col="green",
         xlab="Days since 01-01-2005", ylab="Deseasoned signal",
         main="Deseasoned LAI"
         )
    acf(na.omit(d.af.deseas))
    pacf(na.omit(d.af.deseas))
}

# TODO: write all plots out to pdf
# TODO: detrend using filter(12), deseason with filter(3)
# TODO: try polynomial deseasoning

# Close graphics output file lock
dev.off()

# SOURCES:
# 1. http://geog.uoregon.edu/GeogR/topics/netCDF-read-ncdf4.html
# 2. https://rpubs.com/ryankelly/ts2
