#
library(fdapace)
data2 = read.csv("leafidx.csv")
#Take alook to see how we can fit a time-series model.
test1 = which(data2[,1] == 26.25) #logitude
plot(1:1140, data2[test1[1],-(1:2)], ylim = c(0, 20), type = 'l', col = "grey")
temp = sapply(2:15, function(i) lines(1:1140, data2[test1[i],-(1:2)], col = "grey"))

months = (1:12)*95
plot(t[months], data2[test1[1],months + 2], ylim = c(0, 20), type = 'l', col = "grey")
temp = sapply(2:15, function(i) lines(t[months], data2[test1[i],months + 2], col = "grey"))

range(data2[test1, -(1:2)])
# PCA
# Note: it took a while - large matrix
pca1 = princomp(data2[,-(1:2)])
print(pca1)
plot(pca1)
screeplot(pca1) #turn out first component explains most of the variance.

pc1 = prcomp(data2[,-(1:2)])
plot(pc1)
screeplot(pc1, npcs = 13)
dim(pc1$x)
fit1 = arima(pc1$x[,1]) 
fit1

# fPCA
library(fdapace)
tmplist = MakeFPCAInputs(ID = rep(1:dim(data2)[1], each = 1140), tVec = rep(1:1140, by = dim(data)[1]), 
                         y = as.vector(t(data2[,-(1:2)])))
fpca1 = FPCA(tmplist$Ly, tmplist$Lt, optns = list(dataType = "Dense"))
fpca1$xiEst
fpca1$phi
summary(fpca1)
