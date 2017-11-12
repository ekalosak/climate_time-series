# Author: Eric
# Detrend raw Leaf Area Index observations

library(TTR) # runSD() - running standard dev

source("load_data.R") #, echo=TRUE)

## Extract a single month each year as an attempt to de-season the data
st.mo = 1 # starting month, 1 is Jan 01, 2005
mo.ixs = seq(from=st.mo, to=length(d.af), by=12)
# TODO: try this single season deseasonalization
# TODO: deseasonalize into 12 potentially covarrying timeserieses

## Standardize variance, detrend, deseason etc. standard time series analysis

# detrend with polynomial model
df = data.frame(obs=d.af, time=t)
degree = 15
fit.lm = lm(obs ~ poly(time, degree, raw=TRUE), df)
d.af.detr.nosd = d.af - predict(fit.lm)
# d.af.detr = d.af-predict(fit.lm)

# standardize variance
rsd = runSD(d.af.detr.nosd, 12*10) # 10 year window for sd calculation
msd = lm(rsd~t)
msd.p = predict.lm(msd, data.frame(t))
d.af.detr = d.af.detr.nosd/msd.p

# try using differences
d.af.diff = diff(d.af)

# moving average to de-season
filter.length = 12 # 12 used because t[1]-t[13] = -365
filt.ma = filter(d.af.detr, rep(1,filter.length)/filter.length)

d.af.deseas = filt.ma
d.af.deseas.resid = d.af.detr - filt.ma

# TODO: detrend using filter(12), deseason with filter(3)
# TODO: try polynomial deseasoning

# SOURCES:
# 1. https://rpubs.com/ryankelly/ts2
