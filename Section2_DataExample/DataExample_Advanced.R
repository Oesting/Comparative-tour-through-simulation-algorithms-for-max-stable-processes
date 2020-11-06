## preparations ######################################################################

## load the data 
load("Section2.RData")

## load additional spatial information for plots
library(rnaturalearth)
library(rgdal)
NL_sp <- ne_countries(scale = "medium", country = "Netherlands")
coast_sp <- ne_coastline(scale = "medium")
lakes_sp <- ne_download(scale="medium", type="lakes", category = "physical")

## auxiliary plot functions
source("AuxiliaryFunctions.R")


## view data and extract summer 14 day block maxima ####################################

## plot stations and inland grid
str(stations)
plot.all.stations(display="names")
plot.all.stations(display="numbers")

## extract/plot summer temperatures 
## NOTE: temperatures are measured in 0.1 degree Celsius
str(temperature)

library(lubridate)
dates <- seq(from=ymd(min(temperature$date)), to=ymd(max(temperature$date)), by="1 day" )
summer <- dates[month(dates)==6 & day(dates) > 4 | month(dates)==7 | month(dates)==8 & day(dates) < 28]
summer.temperature <- subset(temperature, date %in% summer)

plot.summer.temperature(stations$stn[1])
plot.summer.temperature(stations$stn[2])

## extract/plot 14 day block maxima of summer temperatures
identical(summer.temperature, summer.temperature[with(summer.temperature, order(stn,date)),])
maxima14days <- summer.temperature[seq(from=1, by=14, length.out = dim(summer.temperature)[1]/14),]
maxima14days$temp <- apply(matrix(summer.temperature$temp,nrow=14),2,max)
maxima14days <- merge(maxima14days,stations,by="stn")

plot.max14days.temperature(stations$stn[1])
plot.max14days.temperature(stations$stn[2])


## estimation of marginal distributions ################################################

## univariate gev fit to 14 day summer maximum temperatures
library(extRemes)
fit.gev <- fevd(x=temp, data=maxima14days[complete.cases(maxima14days),], 
                type="GEV", method = "MLE", location.fun = ~lon + lat + alt)

## estimated gev parameters
mu0 <- fit.gev$results$par[1]
mu1 <- fit.gev$results$par[2]
mu2 <- fit.gev$results$par[3]
mu3 <- fit.gev$results$par[4]
scale <- fit.gev$results$par[5]
shape <- fit.gev$results$par[6]

## sanity check: QQ plots
location.scale.transformed.temp <- with(maxima14days, 
                                        (temp - (mu0 + mu1*lon + mu2*lat + mu3*alt))/scale) 
location.scale.transformed.empirical.quantiles <- sort(location.scale.transformed.temp)
n <- length(location.scale.transformed.empirical.quantiles)
theoretical.quantiles.shape <- qevd((1:n-0.5)/n, loc=0, scale = 1, shape = shape, type="GEV")
plot(theoretical.quantiles.shape, location.scale.transformed.empirical.quantiles,
     xlab="Model quantiles", ylab="Empirical quantiles of transformed data")
abline(0,1)

QQ.plot(stations$stn[1])
QQ.plot(stations$stn[2])

## evaluate estimated GEV location function on the inland grid
location.param <- with(inland.grid, mu0 + mu1*lon + mu2*lat + mu3*elevation_geonames)
inland.grid <- cbind(inland.grid,location.param)
library(viridis)
library(latticeExtra)
## REMINDER: temperatures have been measured in 0.1 degree Celsius
levelplot(location.param/10 ~ lon * lat, inland.grid, 
          panel=panel.levelplot.raster,  col.regions = magma, cuts=50)


## estimation of spatial dependence structure ########################################

## M estimator (for isotropic BR process)
library(tailDepFun)
## ?EstimationBR

## bring temperature data to the right matrix format
identical(maxima14days, maxima14days[with(maxima14days, order(stn,date)),])
maxima.matrix <- matrix(maxima14days$temp,ncol=dim(stations)[1])
maxima.matrix <- maxima.matrix[complete.cases(maxima.matrix),]

## euclidean coordinates of locations
locations <- as.matrix(stations[,c("lon","lat")])
locations[,"lat"] <- lat.factor * locations[,"lat"]
colnames(locations)[2] <- "modif.lat"

## check if euclidean approximation is OK
library(geosphere)
geosphere.distances <- distm(stations[,c("lon","lat")],stations[,c("lon","lat")], fun=distVincentyEllipsoid)
euclidean.distances <- outer(1:18,1:18, function(i,j){sqrt((locations[i,"lon"]-locations[j,"lon"])^2 + (locations[i,"modif.lat"]-locations[j,"modif.lat"])^2)})
plot(geosphere.distances/1e3,euclidean.distances*dist.factor*1e2)
abline(0,1)

## selection of station pairs
plot.pairs(station.pairs)

## execution of the following estimation lasted several hours (i7 CORE, 8th Gen)
if (FALSE) {
  BR.param <- EstimationBR(maxima.matrix, locations, station.pairs, k=50, 
                           method = "Mestimator", isotropic = TRUE, covMat = FALSE, iterate = T)
}
## the result has been saved
BR.param

## sanity check: 

## theoretical bivariate extremal coefficient function
theoretical.ecf <- function(h){2*pnorm(sqrt(2*(h/BR.param$theta[2])^BR.param$theta[1])/2)}

## non-parametric estimates of bivariate extremal coefficients
all.pairs <- matrix(NA,nr=0,nc=dim(stations)[1])
for (i in 2:dim(stations)[1]){
  for (j in 1:(i-1)){
    all.pairs <- rbind(all.pairs, as.integer(stations$stn %in% stations$stn[c(i,j)]))
  }
}
dim(all.pairs)[1]==choose(dim(stations)[1],2)

ec.estimates <- pair.distances <- rep(NA,dim(all.pairs)[1])
col.vector <- rep("black",dim(all.pairs)[1])
library(evd)
for (i in 1:dim(all.pairs)[1]){
  tmp.pair <- which(all.pairs[i,]==1)
  pair.distances[i] <- geosphere.distances[tmp.pair[1],tmp.pair[2]]
  Y=cbind(subset(maxima14days, stn==stations[tmp.pair[1],"stn"])$temp,
          subset(maxima14days, stn==stations[tmp.pair[2],"stn"])$temp)
  ec.estimates[i]  <- 2*abvnonpar(data = Y, method="cfg", empar=T)
  ##
  a <- which(station.pairs[,tmp.pair[1]]==1)
  b <- which(station.pairs[,tmp.pair[2]]==1)
  if(any(a %in% b)){col.vector[i] <-"red"}
}

## comparison of non-parametric estimates of extremal coefficients with model ecf
plot(pair.distances/1e3,ec.estimates, ylim=range(ec.estimates),pch=20,
     ylab="extremal coefficients (CFG-empmar estimates)",
     xlab="distance [km]",col=col.vector)
lines((1:400)/100*dist.factor*1e2,theoretical.ecf((1:400)/100))


## simulation from the fitted max-stable process ###################################

## make simulation algorithms available
source("Algorithms.R") 

## load a fast random number generator
library(dqrng)
rnorm <- dqrnorm

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

## simple marginal check on a sample of 30 locations
sample.grid.locations <- sample(1:dim(simu)[2],30,rep=F)
boxplot(log(simu[,sample.grid.locations])) 
abline(h=qgumbel(c(0.25,0.5,0.75)),col="red")

## transform the simulations to the estimated marginal distributions
no.simu <- dim(simu)[1]
simu.transformed <- simu
for (i in 1:no.simu){
  simu.transformed[i,] <- qgev(pevd(simu[i,], loc=1, scale=1, shape=1,type="GEV"),
                               loc=location.param, scale=scale, shape=shape)
}

## plot of six simulations
no.plot <- 6 ## set this to 1 for one simulation
sample.simulations <- sample(1:dim(simu.transformed)[1], no.plot, rep=F)
data <- data.frame(lon = rep(inland.grid$lon,no.plot),
                   lat = rep(inland.grid$lat,no.plot),
                   simu.tr = as.vector(t(simu.transformed[sample.simulations,]))/10,
                   no.simu = rep(1:no.plot,each=dim(inland.grid)[1]))
levelplot(simu.tr ~ lon * lat | factor(no.simu), data, 
          panel = panel.levelplot.raster, col.regions = magma, cuts=50)


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

## we consider three sub-regions
plot.all.stations()
with(subset(inland.grid,region=="SW"), points(lon,lat,col="yellow",pch="*"))
text(4.8,51.8,labels = "SW",cex=1.5)
with(subset(inland.grid,region=="SE"), points(lon,lat,col="red",pch="*"))
text(5.6,51.8,labels = "SE",cex=1.5)
with(subset(inland.grid,region=="NE"), points(lon,lat,col="green",pch="*"))
text(6.5,52.8,labels = "NE",cex=1.5)

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

## illustration of the joint distribution of the simulated region maxima
library("rgl")
limits <- range(area.maxima)
with(area.maxima,
     rgl.plot3d(SW,SE,NE,"SW","SE","NE",limits,limits,limits,thresh,thresh,thresh,
                col.vector=c("black","red")[1+joint.exceedance], 
                col.vector.proj = c("gray","red")[1+joint.exceedance]))



