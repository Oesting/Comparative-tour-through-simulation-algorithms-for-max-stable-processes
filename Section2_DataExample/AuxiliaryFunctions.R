## Auxiliary plot functions for Section 2

plot.aggregation.2groups <- function(data, VAR, FUN1, VAR1, FUN2, VAR2, aggrFUN=mean, ... ){
  require(fields)
  stopifnot(c(VAR,VAR1,VAR2) %in% colnames(data))
  new.data <- aggregate(data[,VAR], 
                        by = list(Fct1 = FUN1(data[,VAR1]), 
                                  Fct2 = FUN2(data[,VAR2])), 
                        FUN = aggrFUN)
  new.data <- new.data[with(new.data, order(Fct2,Fct1)),]
  one.values <- unique(as.vector(new.data[,"Fct1"]))
  two.values <- unique(as.vector(new.data[,"Fct2"]))
  new.matrix <- matrix(new.data[,"x"],nr=length(two.values),nc=length(one.values),byrow=T)
  two.labels <- as.character(two.values)
  one.labels <- as.character(one.values)
  ltwo <- length(two.labels)
  lone <- length(one.labels)
  fields::image.plot(1:lone,1:ltwo,t(new.matrix),xaxt="n",yaxt="n",xlab="",ylab="",...)
  axis(LEFT<-2, at=1:ltwo, labels=two.labels, las=HORIZONTAL <- 1,cex.axis=0.7)
  axis(BELOW<-1, at=1:lone, labels=one.labels, las=VERTICAL <-2, cex.axis=0.7)
}

plot.all.stations <- function(display="",title.text="Stations"){
  stopifnot(exists("NL_sp") & exists("lakes_sp") & exists("coast_sp"))
  stopifnot(exists("inland.grid") & exists("stations"))
  plot(NL_sp, col="lightgray", 
       xlim=range(inland.grid$lon)+c(-0.1,0.1), ylim=range(inland.grid$lat)+c(-0.1,0.3)) 
  plot(lakes_sp, col="blue", add=TRUE)
  plot(coast_sp, add = TRUE)
  box()
  title(title.text)
  points(inland.grid[,c("lon","lat")],pch=".")
  points(stations[,c("lon","lat")],pch=20)
  if (display=="names") {text(stations[,c("lon","lat")],stations[,"name"],pos=4,cex=0.7,offset = 0.1)}
  if (display=="numbers") {text(stations[,c("lon","lat")],rownames(stations),pos=4,cex=0.7,offset = 0.1)}
}

plot.one.station <- function(station){
  stopifnot(exists("stations"))
  plot.all.stations(title.text = subset(stations, stn==station)[,"name"])
  points(subset(stations, stn==station)[,c("lon","lat")],pch=17,cex=1.5,col="red")
}

plot.summer.temperature <- function(station){
  require(lubridate)
  stopifnot(exists("stations") & exists("summer.temperature"))
  par(mfrow=c(1,2))
  plot.one.station(station)
  plot.aggregation.2groups(data=subset(summer.temperature, stn==station), 
                           VAR="temp", 
                           FUN1=year, VAR1="date", 
                           FUN2=function(d){yday(d)-leap_year(d)-155}, VAR2="date")
  par(mfrow=c(1,1))
}

plot.max14days.temperature <- function(station){
  require(lubridate)
  stopifnot(exists("stations") & exists("maxima14days"))
  par(mfrow=c(1,2))
  plot.one.station(station)
  plot.aggregation.2groups(data=subset(maxima14days, stn==station), 
                           VAR="temp", 
                           FUN1=year, VAR1="date", 
                           FUN2=function(d){(yday(d)-leap_year(d)-156)/14+1}, VAR2="date")
  par(mfrow=c(1,1))
}

QQ.plot <- function(station){
  require(extRemes)
  stopifnot(exists("stations") & exists("maxima14days"))
  stopifnot(exists("location.scale.transformed.temp") & exists("shape"))
  par(mfrow=c(1,2))
  plot.one.station(station)
  x <- sort(location.scale.transformed.temp[maxima14days$stn==station])
  l <- length(x)
  theoretical.quantiles <- qevd((1:l-0.5)/l,loc=0,scale = 1, shape = shape, type="GEV")
  plot(theoretical.quantiles,x,pch=20,
       xlab="Model quantiles",ylab="Empirical quantiles of transformed data")
  abline(0,1)
  s <- matrix(revd(100*l, loc=0,scale = 1, shape = shape, type="GEV"), nc=100)
  points(rep(theoretical.quantiles,100),apply(s,2,sort),col="blue",pch=".")
  points(theoretical.quantiles,x,pch=20,col="red")
  par(mfrow=c(1,1))
}

plot.pairs <- function(stn.pairs){
  stopifnot(exists("stations"))
  stopifnot(dim(stations)[1]==dim(stn.pairs)[2])
  stopifnot(all(sort(unique(as.vector(stn.pairs))) == c(0,1)))
  stopifnot(all(apply(stn.pairs, 1, sum)==2))
  plot.all.stations(display="numbers",title.text="Station pairs")
  for (i in 1:dim(stn.pairs)[1]){
    pair <- which(stn.pairs[i,]==1)
    arrows(x0 = stations[pair[1],"lon"], y0 = stations[pair[1],"lat"], 
           x1 = stations[pair[2],"lon"], y1 = stations[pair[2],"lat"],
           length = 0.01, lwd=0.5)
  }
  points(stations[,c("lon","lat")],col="red",pch=20)
}

rgl_init <- function(new.device = FALSE, bg = "white", width = 840) { 
  require("rgl")
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 30, zoom = 0.9)
}

rgl_add_axes <- function(x, y, z, xlim, ylim, zlim, 
                         xthresh, ythresh, zthresh,
                         thresh.col = "darkgray",
                         axis.col = "black",
                         xlab = "", ylab="", zlab="", 
                         plane.col = "lightgrey") 
{ require("rgl")
  
  rgl.lines(xlim, rep(ylim[1],2), rep(zlim[1],2), color = axis.col, lwd=2)
  rgl.lines(rep(xlim[1],2), ylim, rep(zlim[1],2), color = axis.col, lwd=2)
  rgl.lines(rep(xlim[1],2), rep(ylim[1],2), zlim, color = axis.col, lwd=2)
  
  rgl.lines(xlim, rep(ythresh,2), rep(zlim[1],2), color = thresh.col, lwd=2)
  rgl.lines(xlim, rep(ylim[1],2), rep(zthresh,2), color = thresh.col, lwd=2)
  
  rgl.lines(rep(xthresh,2), ylim, rep(zlim[1],2), color = thresh.col, lwd=2)
  rgl.lines(rep(xlim[1],2), ylim, rep(zthresh,2), color = thresh.col, lwd=2)
  
  rgl.lines(rep(xthresh,2), rep(ylim[1],2), zlim, color = thresh.col, lwd=2)
  rgl.lines(rep(xlim[1],2), rep(ythresh,2), zlim, color = thresh.col, lwd=2)
  
  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(xlim[2], ylim[1], zlim[1]), c(xlim[1], ylim[2], zlim[1]), c(xlim[1], ylim[1], zlim[2]))
  rgl.points(axes, color = axis.col, size = 3)
  
  # Add axis labels
  rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
            adj = c(1, -1), size = 2.5)
  
  # Add plane
  rgl.quads( x = rep(xlim, each = 2), y = rep(ylim[1],4),
             z = c(zlim[1], zlim[2], zlim[2], zlim[1]), color = plane.col, alpha=0.2)
  rgl.quads( x = rep(xlim, each = 2), z = rep(zlim[1],4),
             y = c(ylim[1], ylim[2], ylim[2], ylim[1]), color = plane.col, alpha=0.2)
  rgl.quads( x = rep(xlim[1], each = 4), z = rep(zlim, each=2),
             y = c(ylim[1], ylim[2], ylim[2], ylim[1]), color = plane.col, alpha=0.2)
  
  # Show tick marks
  axis.lab.col <- "black"
  axis3d('x', pos=c( NA, ylim[1], zlim[1] ), col = axis.lab.col)
  axis3d('y', pos=c( xlim[1], NA, zlim[1] ), col = axis.lab.col)
  axis3d('z', pos=c( xlim[1], ylim[1], NA ), col = axis.lab.col)
}

rgl.plot3d <- function(x1,x2,x3,
                       x1lab,x2lab,x3lab,
                       x1lim,x2lim,x3lim,
                       x1thresh, x2thresh, x3thresh, 
                       col.vector, 
                       col.vector.proj){
  rgl_init()
  rgl.spheres(x1, x3lim[1], x2, r = 0.05, color = col.vector.proj)
  rgl.spheres(x1, x3, x2lim[1], r = 0.05, color = col.vector.proj)
  rgl.spheres(x1lim[1], x3, x2, r = 0.05, color = col.vector.proj)
  rgl.spheres(x1, x3, x2, r = 0.1, color = col.vector)
  rgl_add_axes(x1, x3, x2, 
               x1lim, x3lim, x2lim, 
               x1thresh, x3thresh, x2thresh, 
               xlab=x1lab, ylab=x3lab, zlab=x2lab)
  aspect3d(1,1,1)
}
