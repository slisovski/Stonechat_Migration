# library(devtools)
# install_github("SWotherspoon/SGAT")
# install_github("SWotherspoon/BAStag")
# install_github("Slisovski/GeoLight", ref = "Update_2.01")

library(GeoLight) # ver. 2.01
library(SGAT)
library(BAStag)
library(rgeos)
library(MASS)
library(geosphere)
library(maptools)
data(wrld_simpl)
Sys.setenv(tz="GMT")

## read the basic data
# add 1 day to the tagged date since Japanese standard time is GMT+9
dat <- read.csv("gps location.csv")

bird <- data.frame(
  ID = dat$GeolocatorID,
  lon.breed =	dat$Lon,
  lat.breed =	dat$Lat,
  tagged = 	as.POSIXct(as.Date(dat$Day2014)+1,tz="GMT")
)

######################## focal individual
ID <- bird$ID[10]
########################

wd <- getwd()

d.lux <-  read.table(paste(wd,"/data/", ID, ".lux", sep=""), header = FALSE, skip = 25, 
                     col.names = c("Date","Time", "Light"), 
                     colClasses = c("character", "character", "numeric"))
d.lux$Date <- as.POSIXct(strptime(paste(d.lux$Date, d.lux$Time),
                                  "%d/%m/%Y %H:%M:%S", tz="GMT"))
d.lux$Light <- log(d.lux$Light)

offset <- 4
lightImage(d.lux, offset=offset, zlim=c(0,12))

lon.calib <- bird$lon.breed[bird$ID==ID]
lat.calib <- bird$lat.breed[bird$ID==ID]

# plotting the twilight curves for the calibration (captured) site (i.e., Hokkaido)
tm <- seq(d.lux[1,1],d.lux[nrow(d.lux),1], by="day")
rise <- rep(c(TRUE,FALSE),length(tm))
c.dat <- data.frame(Twilight=twilight(rep(tm,each=2),
                                      lon=lon.calib,lat=lat.calib,
                                      rise=rise,zenith=96),Rise=rise)
tsimagePoints(c.dat$Twilight,offset=offset,pch=16,cex=0.1,
              col=ifelse(c.dat$Rise,"dodgerblue","firebrick"))        

# tracking period
# set the starting time at "the top of August", and the end time at the last date of the lux data
d.track.tm <- c(as.POSIXct("2014-08-01",tz="GMT"),as.POSIXct(d.lux$Date[dim(d.lux)[1]],tz="GMT"))
abline(v=d.track.tm,col="green")

# calibration period
# based on our field experiences, we can assume that stonechats stayed at the breeding sites until mid-August
d.calib.tm <- as.POSIXct(c(as.character(bird$tagged[bird$ID==ID]),"2014-08-15"),tz="GMT")
d.calib <- subset(d.lux,(Date>=d.calib.tm[1] & Date<=d.calib.tm[2]))

threshold <- 0.5
# twl <- preprocessLight(d.lux,threshold,offset=offset,lmax=12)
# write.csv(twl,paste(wd,"/results/",ID,"_twl.csv",sep=""),row.names=FALSE)

twl <- read.csv(paste(wd, "/results/", ID, "_twl.csv", sep=""))
twl[,1] <- as.POSIXct(twl[,1],tz="GMT")
twl <- twilightAdjust(twl,300)
  twl <- twl[twl$Twilight<=as.POSIXct("2015-02-17"),]

lightImage(d.lux,offset=offset,zlim=c(0,12))
tsimagePoints(twl$Twilight,offset=offset,pch=16,cex=0.6,
              col=ifelse(twl$Rise,"dodgerblue","firebrick"))


######################## Calibration

twl_calib <- subset(twl,(Twilight>=d.calib.tm[1] & Twilight<=d.calib.tm[2]))

sun <- solar(twl_calib[,1])
z <- refracted(zenith(sun,lon.calib,lat.calib))

twl_t <- twilight(twl_calib[,1],lon.calib,lat.calib,rise=twl_calib[,2],zenith=max(z)+0.1)
twl_dev <- ifelse(twl_calib$Rise, as.numeric(difftime(twl_calib[,1],twl_t,units="mins")),
                  as.numeric(difftime(twl_t,twl_calib[,1],units="mins")))

hist(twl_dev,freq=F,ylim=c(0,0.065),
     main=paste("Stonechat",ID,"Twilight error"),
     xlab="Twilight error (min)")
seq <- seq(0,80,length=100)
fitml_ng <- fitdistr(twl_dev,"log-Normal")
lines(seq,dlnorm(seq,fitml_ng$estimate[1],fitml_ng$estimate[2]),
      col="firebrick",lwd=3,lty=2)

alpha = c(fitml_ng$estimate[1],fitml_ng$estimate[2])

zenith0 <- median(z)
tol <- 0.2
path  <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith0, tol = tol)
path0 <- thresholdLocation(twl$Twilight, twl$Rise, zenith = zenith0, tol = tol)$x
  path$x[is.nan(path0[,2]),2] <- lat.calib

opar <- par(mfrow=c(2,1),mar=c(2,4,1,1)+0.1)
plot(path$time,path$x[,1],type="b",pch=16,cex=0.5,ylab="Lon",xlab='')
abline(h=lon.calib)
plot(path$time,path$x[,2],type="b",pch=16,cex=0.5,ylab="Lat",xlab='')
abline(h=lat.calib)
par(opar)

######################## Initial departure date

tFirst  <- twl[twl$Deleted==F,1][-nrow(twl[twl$Deleted==F,])]
tSecond <- twl[twl$Deleted==F,1][-1]
type=ifelse(twl[twl$Deleted==F,2],1,2)[-nrow(twl[twl$Deleted==F,])]

twl_gl <- data.frame(tFirst=tFirst,tSecond=tSecond,type=type)
cL <- changeLight(twl_gl,quantile=0.9, days=2, summary=F)

siteMap(path$x, site = cL$site, xlim = c(90,150), ylim = c(-10,70))

# In this case site site a-k are deployment site
departure <- max(tSecond[cL$site==which(letters[]=="k")])

## departure time is "2014-10-14"


######################## Initial path

path <- path$x[twl$Twilight>=(departure-4*24*60*60),]
twl  <- subset(twl,Twilight>=(departure-4*24*60*60))
write.csv(path, paste(wd,"/results/",ID,"_thresholdPath.csv",sep=""),row.names=FALSE)

x0 <- path

x0[1:4,1] <- lon.calib
x0[1:4,2] <- lat.calib

opar <- par(mfrow=c(1,1),mar=c(2,4,1,1)+0.1)
plot(x0,type="n")
plot(wrld_simpl,add=T,col="grey95")
box()
lines(x0,col="blue")
points(x0,pch=16,cex=0.5,col="blue")
par(opar)

tFirst  <- twl[twl$Deleted==F,1][-nrow(twl[twl$Deleted==F,])]
tSecond <- twl[twl$Deleted==F,1][-1]
type=ifelse(twl[twl$Deleted==F,2],1,2)[-nrow(twl[twl$Deleted==F,])]

twl_gl <- data.frame(tFirst=tFirst,tSecond=tSecond,type=type)

fixed <- matrix(FALSE,ncol=2,nrow=nrow(twl_gl))
fixed[1:4,] <- cbind(TRUE,TRUE)

cL <- changeLight(twl_gl,quantile=0.75,days=1.5,fixed=fixed,summary=F)

sites <- cL$site

######################## Spatial mask

xlim = c(90,150)
ylim = c(-10,70)

# empty raster
r <- raster(extent(xlim[1],xlim[2],ylim[1],ylim[2]),resolution=0.15)

wrld <- rasterize(wrld_simpl,r)
wrld[] <- ifelse(wrld[]>1,1,NA)

# distance from the coastline (only for the ocean)
d.mask <- distance(wrld)
# plot(d.mask)

# change probability-distance to coastline relationship
sX <- seq(0,max(d.mask[],na.rm=T),by=200)
sY <- 1 + 5*exp(-((sX)/200000)^1)
# plot(sX,sY)

prob <- approxfun(x=sX,y=sY)

prob.mask <- d.mask
prob.mask[] <- prob(d.mask[])

lookup <- function(r,xlim,ylim){
  r <- as.matrix(r)[nrow(r):1,]
  
  xbin <- seq(xlim[1],xlim[2],length=ncol(r)+1)
  ybin <- seq(ylim[1],ylim[2],length=nrow(r)+1)
  
  function(p){
    r[cbind(.bincode(p[,2],ybin),.bincode(p[,1],xbin))]
  }
}

prob <- lookup(prob.mask,xlim=xlim,ylim=ylim)

log.prior <- function(p){
  f <- prob(p)
  ifelse(is.na(f),-1000,log(f))
}


######################## Model paramters

beta0 <- matrix(c(1,0.2,9,0.25),2,2,byrow=T)
# matplot(0:80,cbind(dgamma(0:80,beta0[1,1],beta0[1,2]),
#                    dgamma(0:80,beta0[2,1],beta0[2,2])),
#        type="l",col=c("red","blue"),lty=1,lwd=2,ylab="")

beta <- matrix(c(ifelse(cL$site==0,beta0[2,1],beta0[1,1]),
                 ifelse(cL$site==0,beta0[2,2],beta0[1,2])),ncol=2)

fixedx <- c(rep(TRUE,4),rep(FALSE,nrow(twl)-4))

z0 <- trackMidpts(x0)


######################## Estelle model

model <- thresholdModel(twl$Twilight,twl$Rise,
                        twilight.model="LogNormal",
                        alpha=alpha,beta=beta,
                        logp.x = log.prior, logp.z = log.prior,
                        x0=x0,z0=z0,
                        zenith=quantile(z,probs=0.975),
                        fixedx=fixedx)
proposal.x <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(x0))
proposal.z <- mvnorm(S=diag(c(0.005,0.005)),n=nlocation(z0))

fit <- estelleMetropolis(model,proposal.x,proposal.z,iters=750,thin=20,chains=2)

plot(wrld_simpl,col="grey90",border="grey10",xlim=range(x0[,1]),ylim=range(x0[,2]))
xm <- locationMean(fit$x)
lines(xm,col=rgb(t(col2rgb("cornflowerblue"))/255,alpha=0.9))
points(xm,pch=16,cex=0.8,col=rgb(t(col2rgb("cornflowerblue"))/255,alpha=0.5))
box()

save(fit,file=paste0(wd,"/results/",ID,"_fit.RData"))


######################## Movement residency analysis

zm <- locationSummary(fit$z,time=fit$model$time,collapse=T)

twl.back <- data.frame(Twilight=twilight(fit$model$time[-length(fit$model$time)],zm$'Lon.50%',zm$'Lat.50%',
                                         fit$model$rise[-length(fit$model$time)],zenith=zenith0,iters=5),
                       Rise=fit$model$rise[-length(fit$model$time)])
lightImage(d.lux,offset=offset,zlim=c(0,12))
tsimagePoints(twl.back$Twilight,offset=offset,pch=16,cex=0.8,
              col=ifelse(twl.back$Rise,"dodgerblue","firebrick"))

twl.gl <- data.frame(tFirst=twl.back[-nrow(twl.back),1],
                     tSecond=twl.back[-1,1],
                     type=ifelse(twl.back[,2],1,2)[-nrow(twl.back)])
fixed <- matrix(FALSE,ncol=2,nrow=nrow(twl.gl))
fixed[1:8,] <- cbind(TRUE,TRUE) # the points falling into the vicinity of the breeding site

# changeLight analysis
cL <- changeLight(twl.gl,quantile=0.9,days=1.5, fixed=fixed,summary=F)
mS <- mergeSites(twl.gl,site=cL$site,distThreshold=350,fixed=fixed,degElevation=90-zenith0)
### considering the measurement errors of our geolocation data: distThreshold =  350km

siteMap(cbind(zm$'Lon.50%',zm$'Lat.50%'),site=mS$site,xlim=xlim,ylim=ylim,type='cross',hull=F)


######################## Migration schedule

zm <- locationSummary(fit$z,time=model$time,collapse=T)
zm$site <- c(mS$site,NA)

write.csv(zm,paste(wd,"/results/",ID,"_summary.csv",sep=""),row.names=FALSE)

## Schedule
tm <- as.POSIXct(apply(cbind(zm[,1], zm[,2]), 1 , mean), tz="GMT", origin="1970-01-01")

site <- zm$site

arr <- tm[which(!is.na(site) & !duplicated(site) & site>0)]
dep <- tm[which(!is.na(site) & !duplicated(site, fromLast = T) & site>0)]

out <- data.frame(Site =  letters[1:length(arr)], Arrival = arr, Departure = dep)
if(!is.na(site[1])) out$Arrival[1] <- NA
if(!is.na(site[length(tm)])) out$Departure[nrow(out)] <- NA

out$Days <- round(apply(out, 1, function(x) as.numeric(difftime(x[3], x[2], units = "days"))),2)


out$Lon <- round(aggregate(zm[site>0, c("Lon.50%", "Lat.50%")], by = list(site[site>0]), median)[,2],1)
out$Lat <- round(aggregate(zm[site>0, c("Lon.50%", "Lat.50%")], by = list(site[site>0]), median)[,3],1)

tmp2      <- data.frame(Site = paste0("mig", 1:(nrow(out)-1)), Arrival = out$Departure[-nrow(out)],
                        Departure = out$Arrival[-1], Days = NA, Lon = NA, Lat = NA)
tmp2$Days <- round(apply(tmp2, 1, function(x) as.numeric(difftime(x[3], x[2], units = "days"))),2)

out <- rbind(out, tmp2)
out <- out[order(out[,3]),]

out$Distance <- NA
for(k in 1:nrow(out)) {
  if(is.na(out$Lat[k])){
    out$Distance[k] <- round(distVincentySphere(c(out$Lon[k-1], out$Lat[k-1]),
                                                c(out$Lon[k+1], out$Lat[k+1]))/1000,1)
  }
}

write.csv(out, paste(wd,"/results/",ID,"_schedule.csv",sep=""), row.names=FALSE)


######################## Map for Appendix

png(paste(wd,"/results/",ID,"_MapAppendix.png",sep=""), res = 100, width = 600, height = 650)
par(bty = "n", mar = c(5,5,0,1))
plot(NA, xlim = xlim, ylim = ylim, xlab = "Longitude", ylab = "Latitude", las = 1, cex.axis = 1.2, cex.lab = 1.5)
sh0 <- MapGen2SL("C:/Users/Simeon/Downloads/WCL_1_5000000.dat")
proj4string(sh0) <- proj4string(wrld_simpl)
sh1 <- gIntersection(sh0, as(extent(xlim[1], xlim[2], ylim[1], ylim[2]), "SpatialPolygons"), byid=T)

bn0 <- MapGen2SL("C:/Users/Simeon/Downloads/InternationalBoundaries.dat")
proj4string(bn0) <- proj4string(wrld_simpl)

plot(sh1, add = T, lwd = 1.4)
plot(bn1, add = T, lwd = 1.4)

Site <- ifelse(is.na(site), 0, site)
cols <- c("firebrick", "goldenrod", "forestgreen", "cornflowerblue", "orange", "darkolivegreen2", "cadetblue2",
          "brown2")

points(zm$`Lon.50%`, zm$`Lat.50%`, pch = 16, col = rgb(0.4,0.4,0.4, alpha = 0.4), cex = 1.2)
siteMap(cbind(zm$'Lon.50%',zm$'Lat.50%'), site=Site, lwd = 3, cex = 1.5,
        type='cross', hull=F, add = T, legend = F, col = cols[1:(length(unique(Site))-1)])
dev.off()
