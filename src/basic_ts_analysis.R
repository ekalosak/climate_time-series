# Author: Eric Kalosa-Kenyon
#
# Opens and plots some leaf area index data

rm(list=ls()) # clean workspace

library(stats) # arima()
library(rpart) # rpart() piecewise flat regression tree

source("detrend.R") # loads data and detrends, stds variance, etc.

# Setup parameters (in load_data.R)
# data_dir = "../data/"
# data_file = "lai_Lmon_CCSM4_rcp45_r2i1p1_200601-210012.nc"
# data_path = paste(data_dir, data_file, sep="")
# attr_name = "lai" # leaf area index

# Plotting setup
# plot.filename = "../pdf/basic_analysis.pdf"
# pdf(file=plot.filename)
par(lwd=2)

## Plot the global LAI at time 0
# extract variables to plot
x = lon # 288 longitudes
y = lat # 192 latitudes
z = d[,,7] # all LAI values across the planet at time 7 (July 2005)

png(filename="../img/LAI_global_t0.png")
leaf.colors = function(x){rev(terrain.colors(x))}
filled.contour(x,y,z, color = leaf.colors, asp = 1,
               plot.axes={
                   axis(1); axis(2);
                   points(
                          floor(lon[loc.af[2]]),
                          floor(lat[loc.af[1]])
                          )
               })
title(xlab="Longitude", ylab="Latitude", main="Global leaf cover, July 2005")
dev.off()

## Plot the data for the time series in a single location in sub-saharan
## Africa (i.e. single v=LAI, single s=African, all T=2000 to 2100)
png(filename="../img/pacf_acf_raw.png")
par(mfrow=c(3,1))
plot(time.mo, d.af,
     xlab="Days since Jan 2005",
     ylab="Tree cover (LAI)",
     main="Sub-Saharan Africa, 2005-2100",
     type="l")
acf(d.af, main="Raw data (ACF)")
pacf(d.af, main="Raw data (PACF)")
dev.off()

## Plot standardize variance, detrend, deseasoned data

# fit piecewise flat curve to data
tree = rpart(obs ~ time, data=df)

# plot detrended data and show piecewise fit
png(filename="../img/detrended_acf_pacf.png")
par(mfrow=c(3,1))
plot(df$time, df$obs, type="l")
lines(df$time, predict(tree), col="purple")
plot(time.mo, predict(detr.poly),
     xlab="Time", ylab="Polynomial signal",
     type="l", col="purple",
     ylim=range(d.af))
title(main=paste("Polynomial trend (poly.degree=", degree, ")", sep=""))
lines(time.mo, d.af)
plot(time.mo, d.af.detr.nosd,
     xlab="Time", ylab="Detrended signal",
     type="l", col="blue")
title(main="Detrended")
dev.off()

# plot standardized variance data
png(filename="../img/standatdized_var_LAI.png")
par(mfrow=c(2,1))
plot(time.mo, rsd, main="Running standard deviation 10 year window",
     xlab="Days since Jan 2005", ylab="Std Dev",
     type="l")
lines(time.mo, sd.p)
plot(time.mo, d.af.detr,
     main="Detrended obsv with standardized StDev",
     xlab="Days since Jan 2005", ylab="f(LAI)",
     type="l", col="blue")
dev.off()

# Check acf and pacf with detrended, variance standardized data
png(filename="../img/detrended_stdzd_acf_pacf.png")
par(mfrow=c(3,1))
plot(time.mo, d.af.detr,
     main="Detrended obsv with standardized StDev",
     xlab="Days since Jan 2005", ylab="f(LAI)",
     type="l", col="blue")
acf(d.af.detr, main="Detrended data (ACF)")
pacf(d.af.detr, main="Detrended data (PACF)")
dev.off()

# Plot differenced data
png(filename="../img/differenced_acf_pacf.png")
par(mfrow=c(3,1))
plot(time.mo[1:length(d.af.diff)], d.af.diff, type="l", col="red",
     main="Differenced", xlab="Days since 01-01-2005", ylab="LAI")
acf(d.af.diff, main="Differenced data (ACF)")
pacf(d.af.diff, main="Differenced data (PACF)")
dev.off()

# Plot the moving average deseasonalization and the spectral justification
png(filename="../img/deseasonalization_spectrum.png")
par(mfrow=c(2,1))

m = ar(d.af.detr)
s = spec.ar(m, main=paste("Spectrum of AR(", m$order, ") ~ detrended",
                          sep=""))
abline(v=0.01402806) # largest spectral peak (1/0.014 = 72)
abline(v=0.08333333) # second spectral peak (note: 1/0.0833 = 12)

plot(time.mo, d.af.detr, type="l",
     xlab="Days since 01-01-2005", ylab="LAI", col="blue")
lines(time.mo, filt.ma, type="l", col="red")
title(main="Moving average seasonality (filter.width=1yr)")
abline(v=time.mo[120])
abline(v=time.mo[120+12]) # show a year width on plot
abline(v=time.mo[361])
abline(v=time.mo[361+72]) # show a 5 year trend width on plot
dev.off()

png(filename="../img/deseasonalization.png")
par(mfrow=c(3,1))
plot(time.mo, d.af.deseas,
     type="l", col="green",
     xlab="Days since Jan, 2005", ylab="Deseasoned signal",
     main="Deseasoned LAI")

acf(na.omit(d.af.deseas))
pacf(na.omit(d.af.deseas))
dev.off()

png(filename="../img/deseasonalization_resid.png")
par(mfrow=c(3,1))
plot(time.mo, d.af.deseas.resid,
     type="l", col="green",
     xlab="Days since Jan, 2005", ylab="Deseasoned residuals",
     main="Deseasoned LAI residuals")

acf(na.omit(d.af.deseas.resid))
pacf(na.omit(d.af.deseas.resid))
dev.off()

png(filename="../img/deseasonalization_resid_difference.png")
par(mfrow=c(3,1))
plot(time.mo[1:length(d.af.deseas.resid.diff)], d.af.deseas.resid.diff,
     type="l", col="red",
     xlab="Days since Jan, 2005", ylab="Differenced residuals",
     main="Deseasoned LAI residuals differenced")

acf(na.omit(d.af.deseas.resid.diff))
pacf(na.omit(d.af.deseas.resid.diff))
dev.off()

# Fit arima(0,1,1) as a first pass to the deseasoned data
mod1 = arima(d.af.deseas, c(0,1,1))

## spectral analysis and sinusoidal fit to deseason
par(mfrow=c(1,1))
par(mfrow=c(2,1))

## Fit time series models
# m = arima(d.af.deseas, c(0,1,1))
m = ar(na.omit(d.af.deseas), order.max=4) # looking at the PACF

## simulate some using the model to spot check
png(filename="../img/ar_sim.png")
s = arima.sim(list(order=c(4,0,0), ar=m$ar), n=length(time.mo))
par(mfrow=c(3,1))
plot(time.mo, d.af.deseas,
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
dev.off()

# TODO: examine ts of just jan v that of july, etc.
# TODO: generate resampled CI from..?

# Close graphics output file lock
graphics.off()

# SOURCES:
# 1. https://rpubs.com/ryankelly/ts2
