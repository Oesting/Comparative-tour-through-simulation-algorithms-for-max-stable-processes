## preparations ######################################################################

load("Section2.RData")
str(maxima14days)         ## NOTE: temperatures are measured in 0.1 degree Celsius

## estimation of marginal distributions ################################################

library(extRemes)         ## univariate gev fit to 14 day summer maximum temperatures
fit.gev <- fevd(x=temp, data=maxima14days[complete.cases(maxima14days),], 
                type="GEV", method = "MLE", location.fun = ~lon + lat + alt)

## estimated gev parameters
mu0 <- fit.gev$results$par[1]
mu1 <- fit.gev$results$par[2]
mu2 <- fit.gev$results$par[3]
mu3 <- fit.gev$results$par[4]
scale <- fit.gev$results$par[5]
shape <- fit.gev$results$par[6]

## evaluate estimated GEV location function on the inland grid
location.param <- with(inland.grid, mu0 + mu1*lon + mu2*lat + mu3*elevation_geonames)

## estimation of spatial dependence structure ########################################

library(tailDepFun)       ## M estimator (for isotropic BR process)
## ?EstimationBR

## bring temperature data to the right matrix format
identical(maxima14days, maxima14days[with(maxima14days, order(stn,date)),])
maxima.matrix <- matrix(maxima14days$temp,ncol=dim(stations)[1])
maxima.matrix <- maxima.matrix[complete.cases(maxima.matrix),]

## euclidean coordinates of locations
locations <- as.matrix(stations[,c("lon","lat")])
locations[,"lat"] <- lat.factor * locations[,"lat"]
colnames(locations)[2] <- "modif.lat"

## selection of station pairs
station.pairs 

## execution of the following estimation lasted several hours (i7 CORE, 8th Gen)
if (FALSE) {
  BR.param <- EstimationBR(maxima.matrix, locations, station.pairs, k=50, 
                           method = "Mestimator", isotropic = TRUE, covMat = FALSE, iterate = T)
}
## the result has been saved
BR.param  

## simulation from the fitted max-stable process ###################################

source("Algorithms.R")   ## make simulation algorithms available

## inland coordinates
inland.grid.locations <- data.frame(lon = inland.grid$lon, modif.lat = lat.factor * inland.grid$lat)
coord <- as.matrix(inland.grid.locations)

## semi-variogram 
vario <- function(h) (sqrt(sum(h^2))/BR.param$theta[2])^BR.param$theta[1]

## the following example code for (just) 3 simulations may take 
## a few minutes / an hour (indicative time with i7 CORE, 8th Gen)
if (FALSE){
  simu <- simu_extrfcts(no.simu=3, coord=coord, vario=vario, type="brownresnick")
} 

## instead of simply setting no.simu=60000
## we ran 60 times 1000 simulations on a remote server and merged the output files
## exemplarily we take a look at the first 600 simulations
load("inlandBRsimu_0001.RData")
simu <- simu$res

## transform the simulations to the estimated marginal distributions
no.simu <- dim(simu)[1]
simu.transformed <- simu
library(evd)
for (i in 1:no.simu){
  simu.transformed[i,] <- qgev(pevd(simu[i,], loc=1, scale=1, shape=1,type="GEV"),
                               loc=location.param, scale=scale, shape=shape)
}

## Question 1 #################################################################
## based on inappropriate 600 samples (instead of 60000)

## distribution of the simulated maximum temperature (over a 14 day period) 
## across the entire inland grid
## REMINDER: temperatures have been measured in 0.1 degree Celsius
inland.maxima <- apply(simu.transformed, 1, max)/10
inland.maxima.upper.endpoint <- max(location.param-scale/shape)/10

## return levels if the data were representative of a stationary regime 
## note: 84/14=6 intervals per year
q1 <- quantile(inland.maxima,1-1/(84/14*1e1)) ## level to be seen on average every 10 years
q2 <- quantile(inland.maxima,1-1/(84/14*1e2)) ## level to be seen on average every 100 years
q3 <- quantile(inland.maxima,1-1/(84/14*1e3)) ## level to be seen on average every 1000 years

hist(inland.maxima, breaks=35, freq=F, col="lightgray",
     xlim=c(min(inland.maxima),inland.maxima.upper.endpoint),
     xlab="temperature",ylab="density",
     main="Maximum inland temperature (14 days)")
lines(density(inland.maxima))

abline(v=inland.maxima.upper.endpoint)
abline(v=c(q1,q2,q3),col=c(rgb(0,1,0),rgb(0,0,1),rgb(1,0,0)))

legend("topleft", title="Return level estimates",
       col=c(rgb(0,1,0),rgb(0,0,1),rgb(1,0,0),"black"), lwd=2, cex=0.9,
       legend=c(paste("    10 years: ",(round(as.numeric(q1),1))),
                paste("  100 years: ",(round(as.numeric(q2),1))),
                paste("1000 years: ",(round(as.numeric(q3),1))),
                paste("Upper endp:",(round(inland.maxima.upper.endpoint,1)))))

## Question 2 #################################################################
## based on inappropriate 600 samples (instead of 60000)

## joint distribution of the simulated maximum temperatures (over a 14 day period) 
## across the three subregions of the inland grid
## REMINDER: temperatures have been measured in 0.1 degree Celsius
SW.max <- apply(simu.transformed[,inland.grid$region=="SW"], 1, max)/10
SE.max <- apply(simu.transformed[,inland.grid$region=="SE"], 1, max)/10
NE.max <- apply(simu.transformed[,inland.grid$region=="NE"], 1, max)/10
area.maxima <- data.frame(SW=SW.max,SE=SE.max,NE=NE.max)

## threshold
thresh <- 38
joint.exceedance <- with(area.maxima, SW>thresh & SE > thresh & NE> thresh)
p <- mean(joint.exceedance)
p

## estimated probability of a summer, during which
## there is an exceedance in all three regions SW, SE, NE 
## simultaneously during at least one of the six 14 day intervals
1-(1-p)^6 
6*p

## REMINDER: CAUTION!
## 600 simulations are far too low for a robust result.
## We refer to Section 2 for further reasons why to be cautious.




