# Author: Cody Carroll
#
# Changepoint and Spectral Analyses

source("load_data.R", echo=TRUE)



#Change point
library(ggplot2)
library(ggfortify)
ts.af=ts(d.af)
#autoplot(ts.af, main = "Sub-Saharan Africa LAI")

library(changepoint)
cpt.af=cpt.mean(ts.af, method = "AMOC", class=T, pen.value = 0.01)
summary(cpt.af)
png(filename="../img/changepoint_LAI.png")
autoplot(cpt.af,main = "Sub-Saharan Africa LAI, Changepoint in Mean at 2026")
dev.off()
idx1=1:257
d.af1=d.af[-idx1]
ts.af1=ts(d.af1)
cpt.af1=cpt.mean(ts.af1, method = "AMOC", class=T, pen.value = 0.05)
summary(cpt.af1)
png(filename="../img/changepoint_LAI_seg.png")
autoplot(cpt.af1,main = "Sub-Saharan Africa LAI, First Segment Removed")
dev.off()

png(filename="../img/RawSpec_LAI.png")
raw.spec <- spec.pgram(d.af, main="Periodogram of Raw TS", 
                       col="purple", detrend = F, taper=0.1)
dev.off()



