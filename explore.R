# Author: Eric Kalosa-Kenyon
#
# Opens and plots some leaf area index data

rm(list=ls()) # clean workspace

library(ncdf4) # load .nc data files
library(stats) # arima()

# Setup parameters
data_dir = "./data/"
data_file = "lai_Lmon_CCSM4_rcp45_r2i1p1_200601-210012.nc"
data_path = paste(data_dir, data_file, sep="")
attr_name = "lai" # leaf area index
plot.filename = "basic_analysis.pdf"

# plotting setup
par(lwd=2)

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
filled.contour(x,y,z, color = leaf.colors, asp = 1,
               plot.axes={
                   axis(1); axis(2);
                   points(
                          floor(lon[loc.af[2]]),
                          floor(lat[loc.af[1]])
                          )
               })
title(xlab="Longitude", ylab="Latitude", main="Leaf area, Jan 2005")

## Extract the data for the time series in a single location in sub-saharan
## Africa (i.e. single v=LAI, single s=African, all T=2000 to 2100)
d.af = d[loc.af[2],loc.af[1],] # d.af is african single loc data from 2000 to 2100
plot(t, d.af, xlab="Time", ylab="Tree cover (LAI)",
      main="Single location sub-Saharan African tree cover over a century",
      type="l")

## Extract a single day each year as an attempt to de-season the data
st.day = 1 # starting day, 1 is Jan 01, 2005
# NOTE: the following doesn't work... need to go throught the "t" vector that
# holds the times of each observation
day.ixs = seq(from=st.day, to=length(d.af), by=365)

## Standardize variance, detrend, deseason etc. standard time series analysis

# Check acf and pacf
par(mfrow=c(2,1))
acf(d.af, main="Raw data (ACF)")
pacf(d.af, main="Raw data (PACF)")

# detrend with polynomial model
df = data.frame(obs=d.af, time=t)
degree = 6
fit.lm = lm(obs ~ poly(time, degree, raw=TRUE), df)
d.af.detr = d.af-predict(fit.lm)

par(mfrow=c(2,1))
plot(t, predict(fit.lm),
     xlab="Time", ylab="Polynomial signal",
     type="l", col="purple",
     ylim=range(d.af))
title(main="Polynomial trend")
lines(t, d.af)

plot(t, d.af.detr,
     xlab="Time", ylab="Detrended signal",
     type="l", col="blue")
title(main="Detrended")
par(mfrow=c(1,1))

# Check acf and pacf with detrended data
par(mfrow=c(2,1))
acf(d.af.detr, main="Detrended data (ACF)")
pacf(d.af.detr, main="Detrended data (PACF)")

# try using differences
d.af.diff = diff(d.af)
par(mfrow=c(2,1))
acf(d.af.diff, main="Differenced data (ACF)")
pacf(d.af.diff, main="Differenced data (PACF)")

# moving average to de-season
filter.length = 12 # 12 used because t[1]-t[13] = -365
filt.ma = filter(d.af.detr, rep(1,filter.length)/filter.length)
par(mfrow=c(1,1))
plot(t, d.af.detr, type="l", xlab="Days since 01-01-2005", ylab="LAI")
lines(t, filt.ma, type="l", col="red")
title(main="Moving average seasonality removal")

#d.af.deseas = d.af.detr - filt.ma
d.af.deseas = filt.ma
par(mfrow=c(3,1))
plot(t, d.af.deseas,
     type="l", col="green",
     xlab="Days since 01-01-2005", ylab="Deseasoned signal",
     main="Deseasoned LAI"
     )
acf(na.omit(d.af.deseas))
pacf(na.omit(d.af.deseas))

# TODO: detrend using filter(12), deseason with filter(3)
# TODO: try polynomial deseasoning
# TODO: generate resampled CI from..?

## Fit time series models
# m = arima(d.af.deseas, c(0,1,1))
m = ar(na.omit(d.af.deseas), order.max=4) # looking at the PACF

## simulate some using the model to spot check
ts = m$ar
s = arima.sim(list(order=c(4,0,0), ar=ts), n=length(t))
par(mfrow=c(2,1))
plot(1:length(s), s,
     main="Simulated from AR(4) ~ deseasoned",
     xlab="Months (1140mo=95yr)", ylab="LAI",
     type="l", col="blue")
plot(t, d.af.deseas,
     main="Deseasoned data",
     xlab="Days (3500d=96yr i.e. 2005-2100)", ylab="LAI",
     type="l", col="green")

## spectral analysis
par(mfrow=c(2,1))
plot(1:length(d.af.deseas), d.af.deseas,
     type="l", main="Deseasoned observations", col="green")
sp = spec.ar(m, main="Spectral density of AR(4)",
             col="orange")
par(mfrow=c(1,1))

# Close graphics output file lock
dev.off()

# SOURCES:
# 1. http://geog.uoregon.edu/GeogR/topics/netCDF-read-ncdf4.html
# 2. https://rpubs.com/ryankelly/ts2
