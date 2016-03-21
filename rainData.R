## rainfall data

setwd("~/Dropbox/1_Projects/3_Paper")

library(POT)
library(copula)
library(ADGofTest)
library(foreach)
library(doMC)
library(lattice)
library(data.table)
library(itertools)
#library(maptools)
library(rgdal)
library(ggmap)
#library(mapproj)
## only one should be set to 2
registerDoMC(cores = 1)                                    # just one CPU
MCcores = 2
##

mapply(source, dir("./!functions", rec = TRUE, full.names = TRUE)) 

probs       <- seq(.9, .995, by = 0.001)

## load data
precip <- read.table("./data/1951_2009_daily_values_homogenized.txt")
meta1  <- read.table("./data/metadata_allstations.list", sep = ";")
meta2  <- precip[1 : 3, -1]
time   <- precip[-(1 : 3), 1]
precip <- precip[-(1 : 3), -1 ]


## determine region
waterschap.rg      <- readOGR("./data/Waterschapsgrenzen_28992/","AU_AdministrativeUnitPolygon")
waterschap.sp      <- spTransform(waterschap.rg, CRS("+proj=longlat"))
ValleiEnVeluwe.sp  <- subset(waterschap.sp, waterschap.sp$name == "Vallei & Veluwe")
ValleiEnVeluwe.dat <- fortify(ValleiEnVeluwe.sp) 
precipCoordinates  <- SpatialPoints(coords = cbind(as.numeric(meta2[3, ]), as.numeric(meta2[2, ])), CRS("+proj=longlat"))
prInVEV            <- which(!is.na((precipCoordinates %over% ValleiEnVeluwe.sp)$inspire_id))
points             <- precipCoordinates[prInVEV]

## plot region
source("~/Dropbox/1_Projects/3_Paper/scripts/ChangeInStamenMap.R")

mapImage <- get_map(location = c(lon = 5, lat = 52.5),
                    color = "color",
                    #source = "stamen",
                    #maptype = "watercolor",
                    #source = "google",
                    maptype = "watercolor",
                    zoom = 7)
mapImage <- get_stamenmap()

baseMapNLWater <- ggmap(mapImage)
overView <- baseMap + geom_path(aes(x = long, y = lat), col = "red", size = 1, data = ValleiEnVeluwe.dat) #20.1

pdf(file = "./figs/SelectedWaterschapOverView.pdf", width = 4, height = 4)
print(overView)
dev.off()

mapImage2 <- get_map(location = c(lon = 5.75, lat = 52.25),
                     color = "color",
                     source = "stamen",
                     maptype = "watercolor",
                     zoom = 9)

focusView <- ggmap(mapImage2) + 
  geom_path(aes(x = long, y = lat), col = "red", size = 1, data = ValleiEnVeluwe.dat) +
  geom_point(aes(x = coords.x1, y = coords.x2), col = "black", shape = 3, data = data.frame(points@coords))
print(focusView)
print(focusView + geom_point(aes(x = coords.x1, y = coords.x2), col = "blue", shape = 4, data = data.frame(points[5]@coords)))

pdf(file = "./figs/SelectedWaterschapFocusView.pdf", width = 4, height = 4)
#print(focusView)
print(focusView + geom_point(aes(x = coords.x1, y = coords.x2), col = "blue", shape = 4, data = data.frame(points[5]@coords)))
dev.off()

## determine season
dates <- as.Date(1 : nrow(precip), origin = as.Date("1951-01-01"))
times <- seasonSelector(dates, "summer")

## select subset of the data
pr     <- precip[times, prInVEV]
#points <- precipCoordinates[prInVEV]

##
## declustering
##
pr2 <- apply(pr, 2, windowWhichMax, sepTime = 1)
#pr3 <- pr2
#pr2 <- pr


## threshold stability plot
# single site
mTS1 <- foreach(tau = iter(probs), .combine = "rbind") %dopar% {
  tmpData <- pr2[, 10]
  fitgpd(tmpData, threshold = quantile(tmpData, tau))$param[2]
}
#plot(probs[-(93:96)], mTS1[-(93:96)], xlab = "", ylab = "")

# multiple sites (independent)
mSIVInd <- foreach(tau = iter(probs), .combine = "rbind") %dopar% {
  tmp       <- foreach(s = 1 : ncol(pr), .combine = "c") %do% {
    tmpData <- pr2[, s]
    fitgpd(tmpData, threshold = quantile(tmpData, tau))$param[2]
  }
  return(mean(tmp))
}

pdf(file = "./figs/TSRainIndependentEstimation.pdf", width = 7, height = 4)
plot(probs[-(93:96)], mSIVInd[-(93:96)], xlab = "", ylab = "")
dev.off()
abline(v = probs[61], col = 2, lty = 2, lwd = 2)

# multiple sites (common)
mTSCV <- foreach(tau = iter(probs), .combine = "rbind") %dopar% {
  tmpXpar <- list()
  tmpXpar$S <- ncol(pr2)
  tmpXpar$T <- nrow(pr2)
  tmpXpar$loc <- setQuantiles(pr2, tau)
  excess <- getExcesses(pr2, tmpXpar$loc)
  fpar <- function(p, xpar) {
    loc = xpar$loc
    scale = p[1]*loc
    shape = matrix(p[2], nrow = tmpXpar$T, ncol = tmpXpar$S)
    list(loc = loc, scale = scale, shape = shape)
  }
  tmp <- fitGpdFlex(excess, tmpXpar, fpar, 2, start = c(0.5, 0))
  return(tmp$estimate[2])
}

plot(probs[-(93:96)], mTSCV[-(93:96)], xlab = "", ylab = "")
abline(v = .95, col = 2, lty = 2, lwd = 2)

thresholdStabilityData <- data.table(probs = rep(probs,2), values = c(mTSCV, mSIVInd), estimation = factor(c(rep("common", 96), rep("independent", 96))))
#thresholdStabilityPlot <- ggplot(thresholdStabilityData, aes(x = probs, y = values, colour = estimation)) + geom_line() + xlim(c(0.9, 0.99)) + ylim(c(0.05, 0.2))

ggplot(thresholdStabilityData, aes(x = probs, y = values, colour = estimation)) + geom_line() +
  xlim(c(0.9, 0.99)) + ylim(c(0.05, 0.2)) +
  xlab("$\\tau$") + ylab("$\\hat{\\xi}_{\\mathrm{avg}}$") +
  theme(legend.position = "none")

pdf(file = "./figs/thresholdStabilityRainData.pdf", width = 7, height = 4)
print(thresholdStabilityPlot + xlim(0.9, 0.99) + 
        ylim(0.05, 0.2) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none"))
dev.off()



## mean excess plot
mEVInd <- foreach(tau = iter(probs), .combine = "rbind") %dopar% {
  threshold <- setQuantiles(pr, tau)
  excess    <- excesses(pr, threshold)
  tmp       <- foreach(s = 1 : ncol(pr), .combine = "c") %do% {
    mean(excess[, s], na.rm = TRUE)
  }
  return(c(mean(threshold[1, ]), mean(tmp)))
}
plot(mEVInd[-(95:96), 1], mEVInd[-(95:96), 2], xlab = "", ylab = "")
abline(v = mEVInd[61, 1], col = 2, lty = 2, lwd = 2)
abline(v = mEVInd[41, 1], col = 4, lty = 2, lwd = 2)

##
## GoF statistics
##
xpar <- list()
xpar$S <- ncol(pr)
xpar$T <- nrow(pr)
xpar$fpar  <- function(p, xpar) {
  loc   <- xpar$loc
  scale <- p[1] * loc
  shape <- matrix(p[2], xpar$T, xpar$S)
  list(loc = loc, scale = scale, shape = shape)
}
xpar$nOfP  <- 2
xpar$Start <- c(0.5, 0)

## we have to add slight noise
isNotZero <- pr2 > 0
randomModification <- runif(length(isNotZero[isNotZero]), -0.005, 0.005) 
pr3 <- pr2
pr3[isNotZero] <- pr3[isNotZero] + randomModification 

gofStatCRain <- lapply(probs, getRealGpdGoFStatistic, data = pr3
    ,xpar = xpar, estimationMethod = "common")

B <- 10
ADCRainSingle <- sapply(gofStatCRain, getAverageStatistic, area = B, type = "AD")
KSCRainSingle <- sapply(gofStatCRain, getAverageStatistic, area = B, type = "KS")
critAD95RainSingle.GC <- sapply(bootStatsRain.GC, getQuantileAverageStatistic, area = B, type = "AD", tau = 0.95)
critKS95RainSingle.GC <- sapply(bootStatsRain.GC, getQuantileAverageStatistic, area = B, type = "KS", tau = 0.95)



plot(probs, ADCRainSingle, type = "l")
lines(probs, critAD95RainSingle.GC, col = 2)
plot(probs, KSCRainSingle, type = "l")
lines(probs, critKS95RainSingle.GC, col = 2)




gofStatIRain <- lapply(probs, getRealGpdGoFStatistic, data = pr3
    ,xpar = xpar, estimationMethod = "independent")
KSIRain <- sapply(gofStatIRain, getAverageStatistic, area = 1 : ncol(pr), type = "KS")
ADIRain <- sapply(gofStatIRain, getAverageStatistic, area = 1 : ncol(pr), type = "AD")
plot(probs, ADIRain, type = "l")
plot(probs, KSIRain, type = "l")

B <- 5 #3
KSISingle <- sapply(gofStatIRain, getAverageStatistic, area = B, type = "KS")
ADISingle <- sapply(gofStatIRain, getAverageStatistic, area = B, type = "AD")
critAD95RainSingle.GI <- sapply(bootStatsRain.GI, getQuantileAverageStatistic, area = B, type = "AD", tau = 0.95)
critKS95RainSingle.GI <- sapply(bootStatsRain.GI, getQuantileAverageStatistic, area = B, type = "KS", tau = 0.95)


plot(probs, ADISingle, type = "l")
lines(probs, critAD95RainSingle.GI, col = 2)

KSSingleData <- data.frame(probs = rep(probs,2), values = c(KSISingle, critKS95RainSingle.GI), level = factor(c(rep("real", 96), rep("crit", 96))))
KSSinglePlot <- ggplot(KSSingleData, aes(x = probs, y = values, colour = level)) + geom_line()
print(KSSinglePlot)
pdf(file = "./figs/KSRainSite1.pdf", width = 7, height = 4)
print(KSSinglePlot + scale_color_manual(values = c("red", "black")) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none"))
dev.off()
plot(probs, KSISingle, type = "l")
lines(probs, critKS95RainSingle.GI, col = 2)

##
## tail dependence measurement
##
level <- 0.90
upperTailDep <- foreach(i = 1 : (ncol(pr2) - 1), .combine = "rbind") %do% {
  tmp <- foreach(j = (i + 1) : ncol(pr2), .combine = "rbind") %do% {
    lU <- determineUpperTailDependence(pr2[,i], pr2[, j], level)
    dist <- spDistsN1(points[i], points[j], longlat = TRUE)
    return(c(lU, dist))
  }
  return(tmp)
}

newRhos <- numeric(nrow(upperTailDep))
for(i in 1 : length(newRhos)) {
  newRhos[i] <- computeRhoBasedOnUpperTailDependence(upperTailDep[i,1], 0.9)
}

theoreticalUpperTail <- numeric(length(newRhos))
for(i in 1 : length(newRhos)) {
  theoreticalUpperTail[i] <- computeUpperTailDependence.normalCopula(newRhos[i], 0.99)
}

upperTailDepPar <- matrix(1, nrow = 21, ncol = 21)
foreach(i = 1 : (ncol(pr2) - 1), .combine = "rbind") %do% {
  foreach(j = (i + 1) : ncol(pr2), .combine = "rbind") %do% {
    lU <- determineUpperTailDependence(pr2[,i], pr2[, j], 0.9)
    upperTailDepPar[i,j] <- lU
    upperTailDepPar[j,i] <- lU
  }
}



specParams <- determineGumbelPar(P2p(upperTailDepPar))
specParams <- determineCorrespondingRho(specParams)
specCopula <- tCopula(param = specParams, dim=21, dispstr="un")

newSamples <- 8


upperTailDep <- data.frame("distance" = upperTailDep[,2], "upper tail dependence" = upperTailDep[,1])

fittedTailDep <- loess(upper.tail.dependence ~ distance, data = upperTailDep) #(stat_smooth uses standard loess)
fittedTailDepPar <- matrix(1, nrow = 21, ncol = 21)
foreach(i = 1 : (ncol(pr2) - 1), .combine = "rbind") %do% {
  foreach(j = (i + 1) : ncol(pr2), .combine = "rbind") %do% {
    dist  <- spDistsN1(points[i], points[j], longlat = TRUE)
    index <- which(fittedTailDep$x == dist) 
    lU    <- fittedTailDep$fitted[index]
    fittedTailDepPar[i,j] <- lU
    fittedTailDepPar[j,i] <- lU
  }
}

# fittedParams <- determineGumbelPar(P2p(fittedTailDepPar))
# fittedParams <- determineCorrespondingRho(fittedParams)
# fittedCopula <- tCopula(param = fittedParams, dim=21, dispstr="un")

fittedParams <- P2p(fittedTailDepPar)
for (i in 1 : length(fittedParams)) {
    fittedParams[i] <- computeRhoBasedOnUpperTailDependence(fittedParams[i], 0.9)
}
fittedCopula <- normalCopula(param = fittedParams, dim=21, dispstr="un")
rho <- 0.8336
theta <- 1.933
    
## theoretical tail dependence 
testDep <- fittedTailDepPar[1,5]
nuT     <- 1
rhoT    <- 0.9
rhoT    <- determineCorrespondingRho(determineGumbelPar(testDep))
2*(1 - pt(sqrt(nuT+1)*sqrt(1-rhoT)/sqrt(1+rhoT), df = nuT+1))


testDep  <- 0.75 #as.vector(fittedTailDepPar)
testDep  <- testDep[!testDep==1]
rhoT     <- determineCorrespondingRho(determineGumbelPar(testDep))
nuT     <- 1
testDep2 <- 2*(1 - pt(sqrt(nuT+1)*sqrt(1-rhoT)/sqrt(1+rhoT), df = nuT+1))
sum((testDep-testDep2)^2)


print(depPlot + geom_line(aes(x = distance, y = sorted), data = test2, color = "red"))


depPlot <- ggplot(upperTailDep, aes(x = distance, y = upper.tail.dependence)) + geom_point() + stat_smooth(se=FALSE)

pdf(file = "./figs/TailDependenceOnDistance.pdf", width = 7, height = 4)
print(depPlot + theme(axis.title.x=element_blank(), axis.title.y=element_blank()))
dev.off()


meanUpperTailDep <- mean(upperTailDep$upper.tail.dependence)
theta <- determineGumbelPar(meanUpperTailDep)
tauKendall <- (theta - 1)/theta
rho   <- determineCorrespondingRho(theta)

##
## spatial bootstrap
##
newSamples <- 1000
bootStatsRain.NC <- addGpdGofBootSamples(bootStatsRain.NC, n = newSamples
    ,lambda = probs, data = pr3, xpar = xpar
    ,bootMethod = "copula", copula = normalCopula(param = rho, dim = xpar$S)
    ,estimationMethod = "common", mc.cores = MCcores)
#bootStatsRain.NI <- addGpdGofBootSamples(bootStatsRain.NI, n = newSamples
#    ,lambda = probs, data = prTmp, xpar = xpar
#    ,bootMethod = "copula", copula = normalCopula(param = rho, dim = xpar$S)
#    ,estimationMethod = "independent", mc.cores = MCcores)
bootStatsRain.GC <- addGpdGofBootSamples(bootStatsRain.GC, n = newSamples
    ,lambda = probs, data = pr3, xpar = xpar
    ,bootMethod = "copula", copula = gumbelCopula(param = theta, dim = xpar$S)
    ,estimationMethod = "common", mc.cores = MCcores)
bootStatsRain.GI <- addGpdGofBootSamples(bootStatsRain.GI, n = newSamples
   ,lambda = probs, data = pr3, xpar = xpar
   ,bootMethod = "copula", copula = gumbelCopula(param = 1, dim = xpar$S)
   ,estimationMethod = "independent", mc.cores = MCcores)
#bootStatsRain.UC <- addGpdGofBootSamples(bootStatsRain.UC, n = newSamples
#    ,lambda = probs, data = pr3, xpar = xpar
#    ,bootMethod = "copula", copula = specCopula
#    ,estimationMethod = "common", mc.cores = MCcores)
bootStatsRain.FC <- addGpdGofBootSamples(bootStatsRain.FC, n = newSamples
    ,lambda = probs, data = pr3, xpar = xpar
    ,bootMethod = "copula", copula = fittedCopula
    ,estimationMethod = "common", mc.cores = MCcores)

save(probs, pr, pr2, pr3, xpar, rho, theta, points
    ,bootStatsRain.GC#, bootStatsRain.GI
    ,bootStatsRain.NC#, bootStatsRain.NI
    #,bootStatsRain.UC, specCopula
    ,bootStatsRain.FC, fittedCopula
    #,critAD95Rain.GI, critKS95Rain.GI
    #,critAD95Rain.NI, critKS95Rain.NI
    #,critAD95Rain.GC, critKS95Rain.GC
    #,critAD95Rain.NC, critKS95Rain.NC 
    ,file = "./data/BootstrapRainData.RData")
load("./data/BootstrapRainData.RData")

## check function that writes nothing if error
#for(i in 1 : 96) {
#  bootStatsRain.NC[[i]] <- bootStatsRain.NC[[i]][-301]
#}

# critAD95Rain.GI <- sapply(bootStatsRain.GI, getQuantileAverageStatistic, area = 1 : 21, type = "AD", tau = 0.95)
# critKS95Rain.GI <- sapply(bootStatsRain.GI, getQuantileAverageStatistic, area = 1 : 21, type = "KS", tau = 0.95)
# critAD95Rain.NI <- sapply(bootStatsRain.NI, getQuantileAverageStatistic, area = 1 : 21, type = "AD", tau = 0.95)
# critKS95Rain.NI <- sapply(bootStatsRain.NI, getQuantileAverageStatistic, area = 1 : 21, type = "KS", tau = 0.95)

#critAD95Rain.FC <- sapply(bootStatsRain.FC, getQuantileAverageStatistic, area = 1 : 21, type = "AD", tau = 0.95)
critKS95Rain.FC <- sapply(bootStatsRain.FC, getQuantileAverageStatistic, area = 1 : 21, type = "KS", tau = 0.95)
#critAD95Rain.UC <- sapply(bootStatsRain.UC, getQuantileAverageStatistic, area = 1 : 21, type = "AD", tau = 0.95)
#critKS95Rain.UC <- sapply(bootStatsRain.UC, getQuantileAverageStatistic, area = 1 : 21, type = "KS", tau = 0.95)
#critAD95Rain.GC <- sapply(bootStatsRain.GC, getQuantileAverageStatistic, area = 1 : 21, type = "AD", tau = 0.95)
critKS95Rain.GC <- sapply(bootStatsRain.GC, getQuantileAverageStatistic, area = 1 : 21, type = "KS", tau = 0.95)
#critAD95Rain.NC <- sapply(bootStatsRain.NC, getQuantileAverageStatistic, area = 1 : 21, type = "AD", tau = 0.95)
critKS95Rain.NC <- sapply(bootStatsRain.NC, getQuantileAverageStatistic, area = 1 : 21, type = "KS", tau = 0.95)

#KsIVInd        <- selectThresholdFactorProbMean("KS", bootStatsIVInd.DW, sites = A)
#sUIVIndQ95S16  <- apply(KsIVInd, 1, fnThreshold, z = crit95IVInd.co)

# plot(probs, ADIRain, type = "l")
# lines(probs, critAD95Rain.GI, col = 2)
# lines(probs, critAD95Rain.NI, col = 4)
# 
# plot(probs, KSIRain, type = "l")
# lines(probs, critKS95Rain.GI, col = 2)
# lines(probs, critKS95Rain.NI, col = 4)

gofStatCRain <- lapply(probs, getRealGpdGoFStatistic, data = pr3
                       ,xpar = xpar, estimationMethod = "common")

KSCRain <- sapply(gofStatCRain, getAverageStatistic, area = 1 : ncol(pr), type = "KS")
ADCRain <- sapply(gofStatCRain, getAverageStatistic, area = 1 : ncol(pr), type = "AD")
# #plot(probs, ADCRain, type = "l")
# #plot(probs, KSCRain, type = "l")
# 
# plot(probs, ADCRain, type = "l")
# lines(probs, critAD95Rain.GC, col = 2)
# lines(probs, critAD95Rain.NC, col = 4)

KSData <- data.table(probs = rep(probs, 4), values = c(KSCRain, critKS95Rain.GC, critKS95Rain.NC, critKS95Rain.FC)#, critKS95Rain.UC)
                    ,level = factor(c(rep("1real", 96), rep("2Gumbel", 96), rep("3normal", 96), rep("4n",96))))#, rep("t2",96))))
pdf(file = "./figs/KSRainCommonEstimation.pdf", width = 7, height = 4)
ggplot(KSData, aes(x = probs, y = values, colour = level)) + geom_line() + ylim(0.03, 0.15) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")
dev.off()

plot(probs, KSCRain, type = "l")
lines(probs, critKS95Rain.GC, col = 2)
lines(probs, critKS95Rain.NC, col = 4)




pdf(file = "./figs/TSRainSite5.pdf", width = 7, height = 4)
plot(probs[-(93:96)], mTS1[-(93:96)], xlab = "", ylab = "")
dev.off()






gofStatCRain <- lapply(probs, getRealGpdGoFStatistic, data = pr3
                       ,xpar = xpar, estimationMethod = "common")
statsCRain <- sapply(gofStatCRain, getAverageStatistic, area = 1 : 11, type = "KS")

copulaRain <- gumbelCopula(param = 3.645, dim = sites)
bootStatsCRain.co <- addGpdGofBootSamples(bootStatsCRain.co, n = 50
     ,lambda = probs, data = chosenPrecip2, xpar = xpar
     ,bootMethod = "copula", copula = copulaRain
     ,estimationMethod = "common", mc.cores = MCcores)

crit95CRain.co <- sapply(bootStatsCRain.co, getQuantileAverageStatistic, area = 1 : 11, type = "KS", tau = 0.95)

pdf(file = "./figs/KSRainCommonEstimation.pdf", width = 7, height = 4)
plot(probs, statsCRain, ylim = c(0.04, 0.15), ylab = "", xlab = "")
lines(probs, crit95CRain.co, col = 2)
dev.off()

gofStatIRain <- lapply(probs, getRealGpdGoFStatistic, data = chosenPrecip2
                       ,xpar = xpar, estimationMethod = "independent")


prSingle <- pr2[, 5]
n <- length(prSingle)

params <- foreach(lambda = iter(probs)) %dopar% {
  u      <- as.numeric(quantile(prSingle, probs = lambda))
  peaks  <- getExceedances(prSingle-u, probs[i]) 
  par    <- fitgpd(peaks, threshold = 0)$param
  stat   <- ks.test(peaks[!is.na(peaks)], pgpd, loc = 0, scale = as.numeric(par[1]), shape = as.numeric(par[2]))$statistic
  return(list(loc = as.numeric(u), scale = as.numeric(par[1]), shape = as.numeric(par[2]), stat = stat))
}

result <- foreach(i = 1 : length(params)) %dopar% {
  result <- foreach(j = 1 : 1000, .combine = "c") %do% {
    values <- runif(n)
    peaks  <- getExceedances(values, probs[i])
    peaks  <- qgpdIgnoreNA(peaks, loc = 0, scale = params[[i]]$scale, shape = params[[i]]$shape, lambda = probs[i])
    stat   <- ks.test(peaks[!is.na(peaks)], pgpd, loc = 0, params[[i]]$scale, shape = params[[i]]$shape)$statistic
    return(stat)
  }
  return(result)
}

newValues <- runif(xpar$n)

# bootStatsCRain.co <- addGpdGofBootSamples(bootStatsCRain.co, n = 200
#     ,lambda = probs, data = chosenPrecip2, xpar = xpar
#     ,bootMethod = "copula", copula = copulaRain
#     ,estimationMethod = "common", mc.cores = MCcores)
bootStatsIRain.co <- addGpdGofBootSamples(bootStatsIRain.co, n = 200
    ,lambda = probs, data = chosenPrecip2, xpar = xpar
    ,bootMethod = "copula", copula = copulaRain
    ,estimationMethod = "independent", mc.cores = MCcores)
# 
# save(gofStatCRain, gofStatIRain
#     ,bootStatsCRain.co, bootStatsIRain.co
#     ,file = "./data/StatisticsRain.RData")


statsIRainSingle <- sapply(gofStatIRain, getAverageStatistic, area = 1, type = "KS")
crit95IRainSingle.co <- sapply(bootStatsIRain.co, getQuantileAverageStatistic, area = 1, type = "KS", tau = 0.95)

pdf(file = "./figs/KSRainSite1.pdf", width = 7, height = 4)
plot(probs, statsIRainSingle, ylim = c(0.00, 0.15), ylab = "", xlab = "")
lines(probs, crit95IRainSingle.co, col = 2)
dev.off()



## considerations for simulation paper

dailyPrecipSum <- apply(choosenPrecip, 1, sum)
d10PrecipSum   <- multipleDaySums(dailyPrecipSum, 10)

sortedPrecip <- foreach(s = 1 : dim(choosenPrecip)[2], .combine = "cbind") %do% {
  sort(choosenPrecip[, s], decreasing = TRUE)
}

sum(sortedPrecip[1 : 10, ]) / max(d10PrecipSum) # 4.23 times larger!!!
## how probable is something 2 times larger?
## depends very much on dependence structure

## Is it possible to use the Marshall and Olkin algorithm for simulating copula
## in the following way
## the V which is sampled independently in each step is sampled as Markov Chain?
## => spatial dependence governed by copula
## => temporal dependence only in the generator?
## => Note in the Gumbel copula this should be easily done because stable
## distribution needed for V
## even easier: (bi- trivariate) copula for evolution of V and copula for
## spatial dependence 

##
## individual thresholds
##
u <- setQuantiles(pr, 0.9)
tmp <- getExcesses(pr, u)

fit <- fitgpd(tmp[, 1], 0)
ks.test(tmp[, 1], pgpd, loc = 0, scale = fit$param[1], shape = fit$param[2], exact = FALSE)

