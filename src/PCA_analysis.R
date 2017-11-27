#
#  Author: Amy Kim
# 

### PCA for Spatial-Temporal Data ###

# NOTE:
# "Temporal pattern" explains the dominant temporal variation of time series in all grids, 
# and it is represented by principal components (PCs, a number of time series) of PCA. 
# prcomp(data)$x[,'PC1'] for the most important PC, PC1.
# 
# "Spatial pattern" explains how strong the PCs depend on some variables,
# and it is represented by the loadings of each principal components. 
# prcomp(data)$rotation[,'PC1'].

setwd("~/OneDrive/STA237A/Projects")
library(mapdata)
library(ggplot2) # for map drawing
library(xts)
library(forecast)

print(load('leafidx.RData'))
str(data2)
leafidx_df = as.data.frame(data2)
lat = leafidx_df$lat
lon = leafidx_df$lon
lon2 = ifelse(lon > 180,lon - 360 ,lon )
leafidx = leafidx_df
leafidx$lat = NULL
leafidx$lon = NULL
times = t
leafidx = t(as.data.frame(leafidx)) # columns are where the values belong, rows are the times

drawing <- function(data, map, lonlim = c(-180,180), latlim = c(-90,90)) {
  major.label.x = c("180", "150W", "120W", "90W", "60W", "30W", "0",
                    "30E", "60E", "90E", "120E", "150E", "180")
  major.breaks.x <- seq(-180,180,by = 30)
  minor.breaks.x <- seq(-180,180,by = 10)

  major.label.y = c("90S","60S","30S","0","30N","60N","90N")
  major.breaks.y <- seq(-90,90,by = 30)
  minor.breaks.y <- seq(-90,90,by = 10)
  panel.expand <- c(0,0)

  drawing <- ggplot() +
    geom_path(aes(x = long, y = lat, group = group), data = map) +
    geom_tile(data = data, aes(x = lon, y = lat, fill = val), alpha = 0.3, height = 2) +
    scale_fill_gradient(low = 'white', high = 'red') +
    scale_x_continuous(breaks = major.breaks.x, minor_breaks = minor.breaks.x, labels = major.label.x,
                       expand = panel.expand,limits = lonlim) +
    scale_y_continuous(breaks = major.breaks.y, minor_breaks = minor.breaks.y, labels = major.label.y,
                       expand = panel.expand, limits = latlim) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = 'black'),
          axis.ticks.length = unit(3,"mm"),
          axis.title = element_text(size = 0),
          legend.key.height = unit(1.5,"cm"))

  return(drawing)
}

map.global <- fortify(map(fill=TRUE, plot=FALSE))
dat <- data.frame(lon = lon2, lat = lat, val = leafidx[1,])
sample_plot <- drawing(dat, map.global, lonlim = c(-180,180), c(-90,90))
ggsave("sample_plot.png", sample_plot,width = 6,height=4,units = "in",dpi = 600)

PCAleafidx = prcomp(leafidx, scale = TRUE)
summaryPCAleafidx = summary(PCAleafidx)
summaryPCAleafidx$importance[,c(1,2)] #PC1 34.7% PC2 6.1%
plot(PCAleafidx, main = "Scree Plot")
#plot(PCAdetrend, main = "Scree Plot")
xtable::xtable(summaryPCAleafidx$importance[,1:5], digits = 3)

# Spatial Pattern

loading.PC1 = data.frame(lon = lon2, lat = lat, val = PCAleafidx$rotation[,'PC1']) # loading vectors (eigen vectors)
loading.PC2 = data.frame(lon = lon2, lat = lat, val = PCAleafidx$rotation[,'PC2'])
loading.PC3 = data.frame(lon = lon2, lat = lat, val = PCAleafidx$rotation[,'PC3']) # loading vectors (eigen vectors)
loading.PC4 = data.frame(lon = lon2, lat = lat, val = PCAleafidx$rotation[,'PC4'])
loading.PC5 = data.frame(lon = lon2, lat = lat, val = PCAleafidx$rotation[,'PC5'])
loading.PC6 = data.frame(lon = lon2, lat = lat, val = PCAleafidx$rotation[,'PC6'])

drawing.loadingPC1 <- drawing(loading.PC1, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC1")
drawing.loadingPC2 <- drawing(loading.PC2, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC2")
drawing.loadingPC3 <- drawing(loading.PC3, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC3")
drawing.loadingPC4 <- drawing(loading.PC4, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC4")
drawing.loadingPC5 <- drawing(loading.PC5, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC5")
drawing.loadingPC6 <- drawing(loading.PC6, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC6")

ggsave("loading_PC1.png",drawing.loadingPC1,width = 6,height=4,units = "in",dpi = 600)
ggsave("loading_PC2.png",drawing.loadingPC2,width = 6,height=4,units = "in",dpi = 600)
ggsave("loading_PC3.png",drawing.loadingPC3,width = 6,height=4,units = "in",dpi = 600)
ggsave("loading_PC4.png",drawing.loadingPC4,width = 6,height=4,units = "in",dpi = 600)
ggsave("loading_PC5.png",drawing.loadingPC5,width = 6,height=4,units = "in",dpi = 600)
ggsave("loading_PC6.png",drawing.loadingPC6,width = 6,height=4,units = "in",dpi = 600)



#The "temporal pattern", the first two PC time series, showing the dominant temporal trends of the data

PC1 <- ts(PCAleafidx$x[,'PC1'], frequency = 12) #obs %*% eigenvectors
PC2 <- ts(PCAleafidx$x[,'PC2'], frequency = 12)
PC3 <- ts(PCAleafidx$x[,'PC3'], frequency = 12) #obs %*% eigenvectors
PC4 <- ts(PCAleafidx$x[,'PC4'], frequency = 12)
PC5 <- ts(PCAleafidx$x[,'PC5'], frequency = 12) #obs %*% eigenvectors
PC6 <- ts(PCAleafidx$x[,'PC6'], frequency = 12)


#png("PC-ts.png",width = 6,height = 4,res = 600,units = "in")
plot(as.xts(PC1), ylim = c(-180, 200), main = "PC") # the black one is PC1, major.format = "%Y-%b", type = 'l',
plot(as.xts(PC1), ylim = c(-180, 200), major.format = "%Y-%b", type = 'l', main = "PC")
lines(as.xts(PC2),col='blue',type="l") # the blue one is PC2
lines(as.xts(PC3),col='red',type="l")
lines(as.xts(PC4),col= "green",type="l")
plot(as.xts(PC1), ylim = c(-180, 200), main = "PC") # the black one is PC1, major.format = "%Y-%b", type = 'l',
lines(as.xts(PC4),col='blue',type="l") # the blue one is PC2
lines(as.xts(PC5),col='red',type="l")
lines(as.xts(PC6),col='green',type="l")
dev.off()

par(mfrow = c(3,2))
acf(PC1)
acf(PC2)
acf(PC3)
acf(PC4)
acf(PC5)
acf(PC6)
#pacf
pacf(PC1)
pacf(PC2)
pacf(PC3)
pacf(PC4)
pacf(PC5)
pacf(PC6)
par(mfrow = c(1,1))

pc1fit = auto.arima(PC1)
pc2fit = auto.arima(PC2)

# We can improve the PCA probably by deseasonalizing the data or removing the annual trend by regression, as in the literature suggested. 

### Detrend ###

plot(times, leafidx[,1])
lines(smooth.spline(times, leafidx[,1], cv = TRUE), col = "red")
lines(smooth.spline(times, leafidx[,1], spar = .4), col = "blue")
lines(ksmooth(times,  leafidx[,1], "normal", bandwidth= 200), col=3)
lines(ksmooth(times, leafidx[,1], "box", bandwidth= 400), col=5)

splfit1 = smooth.spline(times, leafidx[,1], cv = TRUE)
res1 = splfit1$yin - splfit1$y


par(mfrow = c(2,1))
plot(times, leafidx[,1],xlab = "Time" , main = "Observations + fitting")
lines(times, splfit1$y, col = "red")
plot(times, res1, type = 'l', xlab= "Time", main = "Detrend")
#lines(times, splfit1$yin, col = "blue")
#plot(times, splfit1$yin, type = "l")
#plot(times, splfit1$y, type = "l")
plot(times, leafidx[,1])
lines(smooth.spline(times, leafidx[,1], cv = TRUE), col = "red")


detrendobs = sapply(1:dim(leafidx)[2], function(i){ fitI <- smooth.spline(times, leafidx[,i], cv = TRUE)
                            resI <- fitI$yin - fitI$y})
#detrendobs[,i] = residual of leafidx[,i]

detrendobs = as.data.frame(detrendobs)

PCAdetrend = prcomp(detrendobs, scale = TRUE)
summaryPCAdetrend = summary(PCAdetrend)
summaryPCAdetrend$importance[,c(1:2)] #PC1 42.6% PC2 9.9%

xtable::xtable(summaryPCAleafidx$importance[,1:6], digits = 3)
xtable::xtable(summaryPCAdetrend$importance[,1:6], digits = 3)

loading.PC1.de = data.frame(lon = lon2, lat = lat, val = PCAdetrend$rotation[,'PC1'])
loading.PC2.de = data.frame(lon = lon2, lat = lat, val = PCAdetrend$rotation[,'PC2'])
loading.PC3.de = data.frame(lon = lon2, lat = lat, val = PCAdetrend$rotation[,'PC3'])
loading.PC4.de = data.frame(lon = lon2, lat = lat, val = PCAdetrend$rotation[,'PC4'])
loading.PC5.de = data.frame(lon = lon2, lat = lat, val = PCAdetrend$rotation[,'PC5'])
loading.PC6.de = data.frame(lon = lon2, lat = lat, val = PCAdetrend$rotation[,'PC6'])

drawing.loadingPC1.de <- drawing(loading.PC1.de, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC1 detrend")
drawing.loadingPC2.de <- drawing(loading.PC2.de, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC2 detrend")
drawing.loadingPC3.de <- drawing(loading.PC3.de, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC3 detrend")
drawing.loadingPC4.de <- drawing(loading.PC4.de, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC4 detrend")
drawing.loadingPC5.de <- drawing(loading.PC5.de, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC5 detrend")
drawing.loadingPC6.de <- drawing(loading.PC6.de, map.global, lonlim = c(-180,180), latlim = c(-90,90)) + ggtitle("PC6 detrend")

ggsave("loading_PC1_de.png",drawing.loadingPC1.de,width = 6,height=4,units = "in",dpi = 600)
ggsave("loading_PC2_de.png",drawing.loadingPC2.de,width = 6,height=4,units = "in",dpi = 600)
ggsave("loading_PC3_de.png",drawing.loadingPC3.de,width = 6,height=4,units = "in",dpi = 600)
ggsave("loading_PC4_de.png",drawing.loadingPC4.de,width = 6,height=4,units = "in",dpi = 600)
ggsave("loading_PC5_de.png",drawing.loadingPC5.de,width = 6,height=4,units = "in",dpi = 600)
ggsave("loading_PC6_de.png",drawing.loadingPC6.de,width = 6,height=4,units = "in",dpi = 600)

PC1.de <- ts(PCAdetrend$x[,'PC1'],frequency = 12)
PC2.de <- ts(PCAdetrend$x[,'PC2'],frequency = 12)
PC3.de <- ts(PCAdetrend$x[,'PC3'],frequency = 12)
PC4.de <- ts(PCAdetrend$x[,'PC4'],frequency = 12)
PC5.de <- ts(PCAdetrend$x[,'PC5'],frequency = 12)
PC6.de <- ts(PCAdetrend$x[,'PC6'],frequency = 12)

png("PC-ts_de.png",width = 6,height = 4,res = 600,units = "in")
par(mfrow = c(1,1))
plot(as.xts(PC1.de), ylim = c(-180, 200), main = "PC Detrend") # the black one is PC1, major.format = "%Y-%b", type = 'l',
lines(as.xts(PC2.de),col='blue',type="l") # the blue one is PC2
lines(as.xts(PC3.de),col='red',type="l")
dev.off()

par(mfrow = c(3,2))
acf(PC1.de)
acf(PC2.de)
acf(PC3.de)
acf(PC4.de)
acf(PC5.de)
acf(PC6.de)
#pacf
pacf(PC1.de)
pacf(PC2.de)
pacf(PC3.de)
pacf(PC4.de)
pacf(PC5.de)
pacf(PC6.de)
par(mfrow = c(1,1))

auto.arima(PC1.de)

pairs(PCAdetrend$x[,1:6])
pairs(PCAleafidx$x[,1:6])

par(mfrow = c(3,2))
ccf(PC1.de, PC2.de)
ccf(PC1.de, PC3.de)
ccf(PC1.de, PC4.de)
ccf(PC2.de, PC3.de)
ccf(PC2.de, PC4.de)
ccf(PC3.de, PC4.de)
par(mfrow = c(1,1))

##Spectral Analysis
pca12.spc = spec.pgram(ts.union(PC1.de, PC2.de), detrend = FALSE)
pca4.spc = spec.pgram(ts.union(PC1.de, PC2.de, PC3.de, PC4.de), detrend = FALSE)
plot(pca4.spc, plot.type = "coherency")
plot(pca4.spc, plot.type = "phase")

# SOURCES: 
# 1. https://stackoverflow.com/questions/41022927/principal-component-analysis-pca-of-time-series-data-spatial-and-temporal-pat

