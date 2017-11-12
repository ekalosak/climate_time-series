# Author: Eric
# Detrend raw Leaf Area Index observations

source("load_data.R", echo=TRUE)

#Setup output file name
plot.filename = "../pdf/detrend_fit_ARIMA.pdf"

# Open outfile in write mode
pdf(file=plot.filename)

## Extract a single day each year as an attempt to de-season the data
st.day = 1 # starting day, 1 is Jan 01, 2005
# NOTE: the following doesn't work... need to go throught the "t" vector that
# holds the times of each observation
day.ixs = seq(from=st.day, to=length(d.af), by=365)

## Standardize variance, detrend, deseason etc. standard time series analysis

# detrend with polynomial model
df = data.frame(obs=d.af, time=t)
degree = 6
fit.lm = lm(obs ~ poly(time, degree, raw=TRUE), df)
d.af.detr = d.af-predict(fit.lm)

# standardize variance
rsd = runSD(d.af.detr.nosd, 12*10) # 10 year window for sd calculation

# try using differences
d.af.diff = diff(d.af)

# moving average to de-season
filter.length = 12 # 12 used because t[1]-t[13] = -365
filt.ma = filter(d.af.detr, rep(1,filter.length)/filter.length)

d.af.deseas = filt.ma
d.af.deseas.resid = d.af.detr - filt.ma

# TODO: detrend using filter(12), deseason with filter(3)
# TODO: try polynomial deseasoning

# Close graphics output file lock
graphics.off()

# SOURCES:
# 1. http://geog.uoregon.edu/GeogR/topics/netCDF-read-ncdf4.html
# 2. https://rpubs.com/ryankelly/ts2
