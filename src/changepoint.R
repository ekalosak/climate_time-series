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
cpt.af=cpt.mean(ts.af)
png(filename="../img/changepoint_LAI.png")
autoplot(cpt.af,main = "Sub-Saharan Africa LAI, Changepoint in Mean at 2026")
dev.off()

raw.spec <- spec.pgram(d.af, main="Periodogram of Raw TS", 
                       col="purple", detrend = F, taper=0.0)


#Fitted spec:
## spectral analysis
png(filename="../img/changepoint_LAI.png")
sp = spec.pgram(na.omit(d.af.deseas), main="Periodogram of Deseasoned TS",
             col="purple", detrend = F)
dev.off()
png(filename="../img/changepoint_LAI.png")
sp = spec.pgram(na.omit(d.af.diff), main="Periodogram of Differenced TS",
                                col="purple", detrend = F) 
dev.off()

