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
df = data.frame(obs=d.af, time=time.mo)
degree = 10
detr.poly = lm(obs ~ poly(time, degree, raw=TRUE), df)
d.af.detr.nosd = d.af - predict(detr.poly)

# standardize variance
# NOTE: the direction of the change in variance may change as the degree of the
# polynomial fit for detrending changes. E.g. try degree=10 and 15
rsd = runSD(d.af.detr.nosd, 12*10) # 10 year window for sd calculation
df$stdv = rsd
# degree2 = 6
# stdv.poly = lm(stdv ~ poly(time, degree2, raw=TRUE), df)
stdv.poly = lm(stdv ~ time, df)
sd.p = predict.lm(stdv.poly, data.frame(time=time.mo))
d.af.detr = d.af.detr.nosd/sd.p
df$detr = d.af.detr

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
