# Author: Cody Carroll
#
# Changepoint and Spectral Analyses

source("explore.R", echo=TRUE)

#Setup parameter
plot.filename = "changepoint_spectral.pdf"


# Open outfile in write mode
pdf(file=plot.filename)


#Change point
library(ggplot2)
library(ggfortify)
ts.af=ts(d.af)
#autoplot(ts.af, main = "Sub-Saharan Africa LAI")

library(changepoint)
cpt.af=cpt.mean(ts.af)
autoplot(cpt.af,main = "Sub-Saharan Africa LAI")
summary(cpt.af)
#cpt.af=cpt.meanvar(ts.af)
#autoplot(cpt.af)

raw.spec <- spec.pgram(d.af, taper = 0)
plot(raw.spec, log = "no", main="No Log")


# Close graphics output file lock
dev.off()