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
?cpt.mean
png(filename="../img/changepoint_LAI.png")
autoplot(cpt.af,main = "Sub-Saharan Africa LAI, Changepoint in Mean at 2026")
dev.off()
cpt.af=cpt.mean(ts.af, method="PELT")
png(filename="../img/multi_changepoint_LAI.png")
autoplot(cpt.af,main = "Sub-Saharan Africa LAI, Changepoint in Mean at 2026")
dev.off()
cpt.af

png(filename="../img/RawSpec_LAI.png")
raw.spec <- spec.pgram(d.af, main="Periodogram of Raw TS", 
                       col="purple", detrend = F, taper=0.1)
dev.off()



