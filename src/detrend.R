# Author: Eric
# Detrend raw Leaf Area Index observations

source("load_data.R", echo=TRUE)

#Setup output file name
plot.filename = "detrend_fit_ARIMA.pdf"

# Open outfile in write mode
pdf(file=plot.filename)

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

d.af.deseas = filt.ma
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

# TODO: detrend using filter(12), deseason with filter(3)
# TODO: try polynomial deseasoning
# TODO: generate resampled CI from..?

# Fit arima(0,1,1) as a first pass to the deseasoned data
mod1 = arima(d.af.deseas, c(0,1,1))

# Close graphics output file lock
dev.off()

# SOURCES:
# 1. http://geog.uoregon.edu/GeogR/topics/netCDF-read-ncdf4.html
# 2. https://rpubs.com/ryankelly/ts2
