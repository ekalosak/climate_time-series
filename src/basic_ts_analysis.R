# Author: Eric Kalosa-Kenyon
#
# Opens and plots some leaf area index data

rm(list=ls()) # clean workspace

library(ncdf4) # load .nc data files
library(stats) # arima()
library(TTR) # runSD() - running standard dev

source("load_data.R")
source("detrend.R")

# Setup parameters
data_dir = "./data/"
data_file = "lai_Lmon_CCSM4_rcp45_r2i1p1_200601-210012.nc"
data_path = paste(data_dir, data_file, sep="")
attr_name = "lai" # leaf area index
plot.filename = "../pdf/basic_analysis.pdf"

# plotting setup
par(lwd=2)

# Open outfile in write mode
pdf(file=plot.filename)

## Plot the global LAI at time 0
# extract variables to plot
x = lon # 288 longitudes
y = lat # 192 latitudes
z = d[,,1] # all LAI values across the planet at time 1

leaf.colors = function(x){rev(terrain.colors(x))}
filled.contour(x,y,z, color = leaf.colors, asp = 1,
               plot.axes={
                   axis(1); axis(2);
                   points(
                          floor(lon[loc.af[2]]),
                          floor(lat[loc.af[1]])
                          )
               })
title(xlab="Longitude", ylab="Latitude", main="LAI at T=1")

plot(t, d.af, xlab="Time", ylab="Tree cover (LAI)",
      main="Single location sub-Saharan African tree cover over a century",
      type="l")

# acf and pacf
par(mfrow=c(2,1))
acf(d.af, main="Raw data (ACF)")
pacf(d.af, main="Raw data (PACF)")

## Plot the data for the time series in a single location in sub-saharan
## Africa (i.e. single v=LAI, single s=African, all T=2000 to 2100)
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
d.af.detr.nosd = d.af-predict(fit.lm)

plot(t, predict(fit.lm),
     xlab="Time", ylab="Polynomial signal",
     type="l", col="purple",
     ylim=range(d.af))
title(main="Polynomial trend")
lines(t, d.af)
plot(t, d.af.detr.nosd,
     xlab="Time", ylab="Detrended signal",
     type="l", col="blue")
title(main="Detrended")

# Check acf and pacf with detrended data
par(mfrow=c(2,1))
acf(d.af.detr, main="Detrended data (ACF)")
pacf(d.af.detr, main="Detrended data (PACF)")

# Plot standardized variance data
plot(t, rsd, main="Running standard deviation 10 year window",
     xlab="Days since Jan 2005", ylab="Std Dev",
     type="l")
msd = lm(rsd~t)
msd.p = predict.lm(msd, data.frame(t))
lines(t, msd.p)
d.af.detr = d.af.detr.nosd/msd.p
plot(t, d.af.detr,
     main="Detrended obsv with stdzd SD",
     xlab="Days since Jan 2005", ylab="f(LAI)",
     type="l")

# Check acf and pacf with detrended data
acf(d.af.detr, main="Detrended data (ACF)")
pacf(d.af.detr, main="Detrended data (PACF)")

# Plot differenced data
plot(t, d.af.detr, type="l", col="blue",
     main="Detrended", xlab="Days since 01-01-2005", ylab="LAI")
plot(t[1:length(d.af.diff)], d.af.diff, type="l", col="red",
     main="Differenced", xlab="Days since 01-01-2005", ylab="LAI")
acf(d.af.diff, main="Differenced data (ACF)")
pacf(d.af.diff, main="Differenced data (PACF)")

# Plot the moving average deseasonalization
par(mfrow=c(1,1))
plot(t, d.af.detr, type="l", xlab="Days since 01-01-2005", ylab="LAI")
lines(t, filt.ma, type="l", col="red")
title(main="One year seasonality removal from detrended observations")

plot(t, d.af.detr, type="l",
     xlab="Days since 01-01-2005", ylab="LAI", col="blue")
lines(t, filt.ma, type="l", col="red")
title(main="Moving average seasonality (filter.width=1yr)")
abline(v=t[120])
abline(v=t[132]) # show a year width on plot

par(mfrow=c(3,1))
plot(t, d.af.deseas,
     type="l", col="green",
     xlab="Days since Jan, 2005", ylab="Deseasoned signal",
     main="Deseasoned LAI")

acf(na.omit(d.af.deseas))
pacf(na.omit(d.af.deseas))

# Fit arima(0,1,1) as a first pass to the deseasoned data
mod1 = arima(d.af.deseas, c(0,1,1))

## spectral analysis and sinusoidal fit to deseason
par(mfrow=c(1,1))
m = ar(d.af.detr)
s = spec.ar(m, main=paste("Spectrum of AR(", m$order, ") ~ detrended",
                          sep=""))
# TODO: find maximal frq -> period
abline(v=0.08) # 0.08 i.e. first peak (note: 1/0.08 = 12.5)
par(mfrow=c(2,1))

## Fit time series models
# m = arima(d.af.deseas, c(0,1,1))
m = ar(na.omit(d.af.deseas), order.max=4) # looking at the PACF

## simulate some using the model to spot check
ts = m$ar
s = arima.sim(list(order=c(4,0,0), ar=ts), n=length(t))
par(mfrow=c(3,1))
plot(t, d.af.deseas,
     main="Deseasoned data",
     xlab="Days (3500d=96yr i.e. 2005-2100)", ylab="LAI",
     type="l", col="green")
plot(1:length(s), s,
     main="Simulated from AR(4) ~ deseasoned",
     xlab="Months (1140mo=95yr)", ylab="LAI",
     type="l", col="blue")
## spectral analysis
sp = spec.ar(m, main="Spectral density of AR(4)",
             col="orange")

# TODO: generate resampled CI from..?

# Close graphics output file lock
graphics.off()

# SOURCES:
# 1. https://rpubs.com/ryankelly/ts2
