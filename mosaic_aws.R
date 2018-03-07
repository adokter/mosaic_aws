#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

##########################################################
#  parsing arguments
##########################################################

if("--help" %in% args){
  cat("mosaic_aws.R [--help] [--suppress-aws-uploads] [--outputdir=<directory>] [--date=<UTC date>]\n")
  cat("   --help                     print this help screen\n")
  cat("   --suppress_uploads         do not store generated imagery in Amazon S3\n")
  cat("   --outputdir=<dir>          optional directory for storing generated imagery\n")
  cat("   --date=<date>              force processing this UTC date\n")
  quit()
}

# determine whether we will upload to S3
S3UPLOAD = !("--suppress_uploads" %in% args)
cat("SCRIPT: ","uploading imagery to S3:",S3UPLOAD,'\n')

# get optional path for storing imagery locally
FIGDIR=args[grep("--outputdir=*",args)]
FIGDIR=sub("--outputdir=","",FIGDIR)
if(length(FIGDIR)!=0) if(!file.exists(FIGDIR)) stop("outputdir not found")

# get forced date if specified
DATEFORCED=args[grep("--date=*",args)]
DATEFORCED=sub("--date=","",DATEFORCED)
if(length(DATEFORCED)!=0){
  d <- try( as.Date(DATEFORCED, format= "%Y-%m-%d %H:%M" ) )
  if( class( d ) == "try-error" || is.na( d ) ) stop("date not in 'YYYY-mm-dd HH:MM' format")
}

# make temporary folder for downloading vp data
VPDIR=paste(tempdir(),"/vp",sep="")
suppressWarnings(dir.create(VPDIR))

##########################################################
#  setup environment
##########################################################

library(bioRad,quietly = T)
library(png,quietly = T)
library(sp,quietly = T)
suppressMessages(library(rgdal,quietly = T))
library(gstat,quietly = T)
suppressWarnings(library(raster,quietly = T))
suppressMessages(library(fields,quietly = T))
library(aws.s3,quietly = T)
suppressMessages(suppressWarnings(library(dplyr,quietly = T)))
library(geosphere,quietly = T)
library(gstat, quietly = T)
if(length(FIGDIR)==0){
  setwd("/opt")
}

Sys.setenv(TZ="UTC")

# S3 bucket where profiles and imagery are stored
BUCKET="vol2bird"

if(file.exists("~/.aws/credentials")){
  getkeys=function() read.csv("~/.aws/credentials",sep="=",stringsAsFactors = F,col.names = "key",strip.white=T)
  Sys.setenv(AWS_ACCESS_KEY_ID=getkeys()["aws_access_key_id",])
  Sys.setenv(AWS_SECRET_ACCESS_KEY=getkeys()["aws_secret_access_key",])
}

# run below lines when testing locally:
# Sys.setenv(MOSAIC_DATE=as.character("2017-10-24 00:40"))
# setwd("~/git/mosaic_aws")

DATE=as.POSIXct(if(Sys.getenv("MOSAIC_DATE")=="") as.character(Sys.time()) else Sys.getenv("MOSAIC_DATE"))
if(length(DATEFORCED)!=0) DATE=as.POSIXct(DATEFORCED)

# UTC hour at which we transition to a new day for grouping radar files
HOUR_NEW_DAY=22

evening.palette <- colorRampPalette(c("#0a0a0a","#1775cb","white"),bias=1.0)(100)

# minimum number of radars for outputting an image
RADARS_MIN=120

##########################################################
#  functions
##########################################################

ipol.vplist = function(vplist,alpha=NA,log=F,alt.min=0, alt.max=Inf,quantity="mtr",method="idw",variogram='auto'){
  if(!(method %in% c("idw","krige"))) stop("unknown interpolation method")
  if(!(quantity %in% c("mtr","u","v","ff","dd","HGHT"))) stop("unknown quantity")
  radar=sapply(vplist, function(x) x$radar)
  locs=radarInfo[sapply(radar,function(x) which(radarInfo$radar==x)),c("lat","lon","height")]
  mtrs=sapply(vplist,mtr)
  # to avoid NA's after taking logarithm
  mtrs[mtrs<1]=1
  # bind the data with radar locations
  data=cbind(radar,locs,mtrs)
  coordinates(data)<-~lon+lat
  proj4string(data)=proj4string(lower48wgs)
  # convert to mercator projection
  data=spTransform(data, CRS("+proj=merc"))
  # interpolate the migration traffic data
  if(method=="idw"){
    if(log) ipol=idw(formula = log10(mtrs) ~ 1, locations = data, newdata = r,nmax=5,idp=2)
    else ipol=idw(formula = mtrs ~ 1, locations = data, newdata = r,nmax=5,idp=2)
  }
  else{
    if(log) g <- gstat(id="var1", formula=log10(mtrs) ~ 1, data=data)
    else g <- gstat(id="var1", formula=mtrs ~ 1, data=data)
    # interpolate the densities
    vario.mtr=variogram(g, map=FALSE)
    if(inherits(variogram,"variogramModel")){
      vario.fit=variogram
    }
    else{
      vario.fit <- fit.variogram(vario.mtr, model=vgm(model='Sph'))
    }
    g <- gstat(g, id="var1", model=vario.fit )
    ipol <- predict(g, newdata=r)
  }
  attributes(ipol)$date=vplist[[1]]$datetime
  attributes(ipol)$radar=vplist[[1]]$radar
  attributes(ipol)$log=log
  attributes(ipol)$data.vp=data
  ipol
}

plotidw <- function(idw,zlim=c(3.5,6), linecol="black",bg="white", file=NA, closedev=T, date.lab='auto',legend.lab="Migration traffic rate [birds/km/h]",col=evening.palette,...){
  if(!is.na(file)) jpeg(file=file, width = 3.5, height = 2.5,units="in",res=300,quality=90)
  if(attributes(idw)$log){
    axislabels=c(0.05,0.5,2,10,50)
    axisticks=log10(1000*axislabels)
    zlim=log10(zlim)
  }
  else{
    axislabels=zlim
    axisticks=axislabels
  }

  bg.start=par()$bg
  fg.start=par()$fg
  par(bg=bg,fg=linecol,col.axis=linecol,col.lab=linecol,col.main=linecol,col.sub=linecol,cex=.1)
  idw@data$var1.pred[idw$var1.pred < zlim[1]]=zlim[1]
  idw@data$var1.pred[idw$var1.pred > zlim[2]]=zlim[2]
  #plot(idw[idw.idx,],add.axis=T,axis.pos=1,scale.shrink=1,zlim=zlim, bg=bg, ...)
  plot(idw[idw.idx,],what="image",zlim=zlim, bg=bg, col=col, ...)
  plot(lower48,add=T,border=linecol,lwd=0.3)
  #plot(cmtdata,add=T,col=linecol,lwd=0.1,pch='.')

  usr <- par("usr")
  xwidth <- usr[2] - usr[1]
  yheight <- usr[4] - usr[3]

  # plot date
  if(date.lab=='auto') date.lab=attributes(idw)$date
  my.date.lab=tryCatch(format(date.lab,"%d %B %Y %H:%M"),error= function(err) {date.lab})
  text(usr[1] + xwidth/23, usr[3] + yheight/5, my.date.lab, col=linecol,cex = 6, font=2,pos=4)

  # plot logo
  if(!is.na(file)){
    rasterImage(logo, usr[2]+usr[2]/5.5, usr[3] + (yheight * 0.1) + yheight/20,
                usr[2] +usr[2]/5.5 - (((usr[3]) - (usr[3] - (yheight * 0.1)))*(res[2]/res[1])),
                usr[3] + yheight/20, interpolate=FALSE,angle=180)
  }

  bins=seq(zlim[1],zlim[2],length.out=length(col)+1)



  # plot radar locations of active radars
  plot(attributes(idw)$data.vp,add=T,col='green',cex=1.5,pch=16)
  # also plot inactive radars
  radarsOffline=radarInfo[!(1:nrow(radarInfo) %in% as.numeric(rownames(attributes(idw)$data.vp@data))),]
  if(nrow(radarsOffline)>0){
    coordinates(radarsOffline)<-~lon+lat
    proj4string(radarsOffline)=proj4string(lower48wgs)
    radarsOffline=spTransform(radarsOffline, CRS("+proj=merc"))
    plot(radarsOffline,add=T,col='red',cex=1.5,pch=16)
  }

  image.plot(idw,
             col=col,
             zlim=zlim,
             #smallplot=c(0.014,0.024,0.12,.82), left aligntment
             smallplot=c(0.05,0.4,0.15,0.165), # bottom allignment
             legend.only=TRUE,
             legend.shrink = 1,
             horizontal=T,
             legend.width=2.5,
             breaks=bins,
             axis.args=list(at=axisticks,fg=linecol, labels=axislabels,
                            col.axis=linecol,
                            cex.axis=4,pos=1,mgp = c(0, 3, 0),
                            lwd.ticks=2),
             legend.args=list(text=legend.lab, side=1, cex=.4, col=linecol, line=8)
  )
  par(bg=bg.start,fg=fg.start,col.axis=fg.start,col.lab=fg.start,col.main=fg.start,col.sub=fg.start)
  if(!is.na(file) && closedev==T) dev.off()
}


plot_wind_barbs = function(cx, cy, direction = 0, speed = NA, arrow=F, fill = rep(0, length(cx)), circle = FALSE, cex = 1, lwd=0.5, col = "black")
{
  ### press is actually height in upper air ###
  ns = length(cx)
  if (length(cy) != ns) stop("X AND Y COORDINATES SHOULD HAVE SAME LENGTH!")

  msg = "ALL VARIABLES SHOULD HAVE SAME LENGTH AS COORDINATES, OR BE MISSING!!!"
  if (ns > 1) {
    if (length(direction) == 1) if (!is.na(direction)) stop(msg)
    if (length(speed) == 1) if (!is.na(speed)) stop(msg)
    if (length(fill) == 1) if (!is.na(fill)) stop(msg)
    if (length(direction) > 1 & length(direction) != ns) stop(msg)
    if (length(speed) > 1 & length(speed) != ns) stop(msg)
    if (length(fill) > 1 & length(fill) != ns) stop(msg)
  }

  if(arrow) direction = direction + 180

  tpar = par()
  size = tpar$csi
  scalex = (tpar$usr[2] - tpar$usr[1]) / tpar$pin[1]
  scaley = (tpar$usr[4] - tpar$usr[3]) / tpar$pin[2]
  scalex = (cex * (scalex * size)) / 5
  scaley = (cex * (scaley * size)) / 5
  for (i in 1 : ns) {
    x = cx[i]
    y = cy[i]
    if (is.na(x) || is.na(y)) next
    spd = speed[i]

    if (circle) {
      ts = seq(0, 2 * pi, length.out = 200)
      RX = sin(ts) * scalex
      X1 = RX + x
      RY = cos(ts) * scaley
      Y1 = RY + y
      if (!is.na(spd)) {
        if (spd == 0) {
          lines(RX * 2 + x, RY * 2 + y, col = col)
        }
      }
      if (fill[i] > 0) {
        lim = c(51, 101, 151, 200)
        polygon(c(x, X1[1 : lim[fill[i]]]), c(y, Y1[1 : lim[fill[i]]]), density = -1, col = col)
      }
      lines(RX + x, RY + y, col = col)
    } #end of circle

    if (!is.na(spd)) {
      if (spd > 0) {
        X1 = 0
        X2 = 0
        Y1 = 0
        Y2 = 5
        if(arrow){
          X1 = c(X1, 0)
          X2 = c(X2, 1, 0, -1)
          Y1 = c(Y1, 5)
          Y2 = c(Y2, 4 ,5, 4)
        }
        else{
          if (spd >= 5 & spd < 10) {
            X1 = c(X1, 0)
            X2 = c(X2, 1)
            Y1 = c(Y1, 5)
            Y2 = c(Y2, 5)
          }
          if (spd >= 10 & spd < 15) {
            X1 = c(X1, 0)
            X2 = c(X2, 2)
            Y1 = c(Y1, 5)
            Y2 = c(Y2, 5)
          }
          if (spd >= 15 & spd < 20) {
            X1 = c(X1, 0, 0)
            X2 = c(X2, 1, 2)
            Y1 = c(Y1, 4, 5)
            Y2 = c(Y2, 4, 5)
          }
          if (spd >= 20 & spd < 25) {
            X1 = c(X1, 0, 0)
            X2 = c(X2, 2, 2)
            Y1 = c(Y1, 4, 5)
            Y2 = c(Y2, 4, 5)
          }
          if (spd >= 25 & spd < 30) {
            X1 = c(X1, 0, 0, 0)
            X2 = c(X2, 1, 2, 2)
            Y1 = c(Y1, 3, 4, 5)
            Y2 = c(Y2, 3, 4, 5)
          }
          if (spd >= 30 & spd < 35) {
            X1 = c(X1, 0, 0, 0)
            X2 = c(X2, 2, 2, 2)
            Y1 = c(Y1, 3, 4, 5)
            Y2 = c(Y2, 3, 4, 5)
          }
          if (spd >= 35 & spd < 40) {
            X1 = c(X1, 0, 0, 0, 0)
            X2 = c(X2, 1, 2, 2, 2)
            Y1 = c(Y1, 2, 3, 4, 5)
            Y2 = c(Y2, 2, 3, 4, 5)
          }
          if (spd >= 40 & spd < 45) {
            X1 = c(X1, 0, 0, 0, 0)
            X2 = c(X2, 2, 2, 2, 2)
            Y1 = c(Y1, 2, 3, 4, 5)
            Y2 = c(Y2, 2, 3, 4, 5)
          }
          if (spd >= 45 & spd < 50) {
            X1 = c(X1, 0, 0, 0, 0, 0)
            X2 = c(X2, 1, 2, 2, 2, 2)
            Y1 = c(Y1, 1, 2, 3, 4, 5)
            Y2 = c(Y2, 1, 2, 3, 4, 5)
          }
          if (spd >= 50 & spd < 55) {
            X1 = c(X1, 0, 0)
            X2 = c(X2, 2, 2)
            Y1 = c(Y1, 4, 5)
            Y2 = c(Y2, 4.5, 4.5)
          }
          if (spd >= 55 & spd < 60) {
            X1 = c(X1, 0, 0, 0)
            X2 = c(X2, 1, 2, 2)
            Y1 = c(Y1, 3, 4, 5)
            Y2 = c(Y2, 3, 4.5, 4.5)
          }
          if (spd >= 60 & spd < 65) {
            X1 = c(X1, 0, 0, 0)
            X2 = c(X2, 2, 2, 2)
            Y1 = c(Y1, 3, 4, 5)
            Y2 = c(Y2, 3, 4.5, 4.5)
          }
        }
        dir = (direction[i] / 360) * 2 * pi
        rot = cbind(c(cos(dir), -sin(dir)), c(sin(dir), cos(dir)))
        S1 = rbind(X1, Y1)
        S2 = rbind(X2, Y2)
        S1 = rot %*% S1
        S2 = rot %*% S2
        S1 = S1 * c(scalex, scaley) + c(x, y)
        S2 = S2 * c(scalex, scaley) + c(x, y)
      }
      if (spd > 0) {
        segments(S1[1, ], S1[2, ], S2[1, ], S2[2, ], col = col, lwd = lwd)
      }
    } #end of (!is.na(spd))
  } #end of ns
  invisible()
}

add_barbs=function(vplist,mtr.min=500,radar.col="white"){
  radar=sapply(vplist, function(x) x$radar)
  locs=radarInfo[sapply(radar,function(x) which(radarInfo$radar==x)),c("lat","lon","height")]
  mtr=sapply(vplist,mtr)
  u=sapply(vplist, function(x) sum(fetch(x,"u")*fetch(x,"dens"),na.rm=T)/sum(fetch(x,"dens"),na.rm=T))
  v=sapply(vplist, function(x) sum(fetch(x,"v")*fetch(x,"dens"),na.rm=T)/sum(fetch(x,"dens"),na.rm=T))
  ff=sqrt(u^2+v^2)
  dd=(0.5*pi-atan2(v,u))*180/pi
  barbdata=cbind(locs,dd,ff,mtr)
  barbdata=barbdata[barbdata$mtr>mtr.min,]
  if(nrow(barbdata)==0) return()
  barbdata[is.na(barbdata)] <- 0
  coordinates(barbdata)<-~lon+lat
  proj4string(barbdata)=proj4string(lower48wgs)
  barbdata=spTransform(barbdata, CRS("+proj=merc"))
  plot_wind_barbs(barbdata@coords[,1],barbdata@coords[,2],180+barbdata$dd,3.5+0*barbdata$ff,col="#fca71d",cex=.35,lwd=0.75,arrow=T)
}

s3file=function(radar,date,bucket=BUCKET,prefix="output",dt=900){
  keyprefix1=paste(prefix,strftime(date,format="/%Y/%m/%d/"),radar,"/",radar,strftime(date+dt,format="%Y%m%d_%H"),sep="")
  keyprefix2=paste(prefix,strftime(date-dt,format="/%Y/%m/%d/"),radar,"/",radar,strftime(date-dt,format="%Y%m%d_%H"),sep="")
  if(keyprefix1==keyprefix2) keys=get_bucket_df(bucket,keyprefix1,max=Inf)
  else keys=rbind(get_bucket_df(bucket,keyprefix1,max=Inf),get_bucket_df(bucket,keyprefix2,max=Inf))
  keys$datetime=as.POSIXct(substr(basename(keys$Key),5,19),format="%Y%m%d_%H%M%S",tz="UTC")
  if(nrow(keys)>0){
    keys=keys[which.min(abs(keys$datetime-date)),]
    if(as.numeric(abs(keys$datetime-date),units="secs")>dt) keys=keys[NULL,]
  }
  keys
}

##########################################################
# calculate day/night terminator
# from https://github.com/JoGall/terminator/blob/master/terminator.R
##########################################################

# ADDS A TERMINATOR TO A MAP TO SHOW DAYTIME / NIGHTTIME REGIONS
# Returns a dataframe of latitude and longitude for the line that separates illuminated day and dark night for any given time
# This is just a port of the Javascript Leaflet.Terminator plugin (https://github.com/joergdietrich/Leaflet.Terminator/blob/master/L.Terminator.js)

rad2deg <- function(rad) {
  (rad * 180) / (pi)
}

deg2rad <- function(deg) {
  (deg * pi) / (180)
}

getJulian <- function(time) {
  # get Julian day (number of days since noon on January 1, 4713 BC; 2440587.5 is number of days between Julian epoch and UNIX epoch)
  (as.integer(time) / 86400) + 2440587.5
}

getGMST <- function(jDay) {
  # calculate Greenwich Mean Sidereal Time
  d <- jDay - 2451545.0
  (18.697374558 + 24.06570982441908 * d) %% 24
}

sunEclipticPosition <- function(jDay) {
  # compute the position of the Sun in ecliptic coordinates
  # days since start of J2000.0
  n <- jDay - 2451545.0
  # mean longitude of the Sun
  L <- 280.460 + 0.9856474 * n
  L = L %% 360
  # mean anomaly of the Sun
  g <- 357.528 + 0.9856003 * n
  g = g %% 360
  # ecliptic longitude of Sun
  lambda <- L + 1.915 * sin(deg2rad(g)) + 0.02 * sin(2 * deg2rad(g))
  # distance from Sun in AU
  R <- 1.00014 - 0.01671 * cos(deg2rad(g)) - 0.0014 * cos(2 * deg2rad(g))

  data.frame(lambda, R)
}

eclipticObliquity <- function(jDay) {
  # compute ecliptic obliquity
  n <- jDay - 2451545.0
  # Julian centuries since J2000.0
  T <- n / 36525
  # compute epsilon
  23.43929111 -
    T * (46.836769 / 3600
         - T * (0.0001831 / 3600
                + T * (0.00200340 / 3600
                       - T * (0.576e-6 / 3600
                              - T * 4.34e-8 / 3600))))
}

sunEquatorialPosition <- function(sunEclLng, eclObliq) {
  # compute the Sun's equatorial position from its ecliptic position
  alpha <- rad2deg(atan(cos(deg2rad(eclObliq)) *
                          tan(deg2rad(sunEclLng))))
  delta <- rad2deg(asin(sin(deg2rad(eclObliq))
                        * sin(deg2rad(sunEclLng))))

  lQuadrant  = floor(sunEclLng / 90) * 90
  raQuadrant = floor(alpha / 90) * 90
  alpha = alpha + (lQuadrant - raQuadrant)

  data.frame(alpha, delta)
}

hourAngle <- function(lng, sunPos, gst) {
  # compute the hour angle of the sun for a longitude on Earth
  lst <- gst + lng / 15
  lst * 15 - sunPos$alpha
}

longitude <- function(ha, sunPos) {
  # for a given hour angle and sun position, compute the latitude of the terminator
  rad2deg(atan(-cos(deg2rad(ha)) / tan(deg2rad(sunPos$delta))))
}

terminator <- function(time, from = -180, to = 180, by = 0.1) {
  # calculate latitude and longitude of terminator within specified range using time (in POSIXct format, e.g. `Sys.time()`)
  jDay = getJulian(time)
  gst = getGMST(jDay)

  sunEclPos = sunEclipticPosition(jDay)
  eclObliq = eclipticObliquity(jDay)
  sunEqPos = sunEquatorialPosition(sunEclPos$lambda, eclObliq)

  lapply(seq(from, to, by), function(i) {
    ha = hourAngle(i, sunEqPos, gst)
    lon = longitude(ha, sunEqPos)
    data.frame(lat = i, lon)
  }) %>%
    plyr::rbind.fill()
}

# # EXAMPLES
# # terminator for current time on world map
#
# # add terminator at current time to world map as shaded region using `geom_ribbon``
# ggplot2::ggplot() +
#   borders("world", colour = "gray90", fill = "gray85") +
#   geom_ribbon(data = terminator(Sys.time(), -180, 190), aes(lat, ymax = lon), ymin = 90, alpha = 0.2) +
#   coord_equal() +
#   ggthemes::theme_map()
#
# # add terminator at a specific time to map of Europe, using a `coord_*` function to crop after drawing shaded region with `geom_ribbon`
# ggplot2::ggplot() +
#   borders("world", colour = "gray90", fill = "gray85") +
#   geom_ribbon(data = terminator(as.POSIXct("2018-01-01 07:00:00 GMT"), -180, 190, 0.1), aes(lat, ymax = lon), ymin = 90, alpha = 0.2) +
#   coord_equal(xlim = c(35, -12), ylim = c(35, 72), expand = 0) +
#   ggthemes::theme_map()

get_terminator=function(date){
  data=terminator(date, -180, 180, 0.1)
  # the terminator function grabbed from https://github.com/JoGall/terminator/blob/master/terminator.R has
  # latitude and longitude erroneously swapped.
  data$lon2=data$lat
  data$lat2=data$lon
  data$lon=data$lon2
  data$lat=data$lat2
  data=data %>% filter(lat >= 29.5 & lat <= 49 & lon > -150 & lon < -49)
  if(!nrow(data)>2) return(NULL)
  coordinates(data)<-~lon+lat
  proj4string(data)=proj4string(lower48wgs)
  # convert to mercator projection
  data=spTransform(data, CRS("+proj=merc"))
  data
}

##########################################################
#  load basemap (code commented out, load from .RData file)
##########################################################
#load("~/Dropbox/radar/NEXRAD/occultation/NEXRAD_ground_antenna_height.RData")
#map=readOGR("./bgmap")
#lower48wgs=map[map$iso_a2=="US" & map$name != "Alaska" & map$name != "Hawaii",]
#lower48=spTransform(lower48wgs, CRS("+proj=merc"))
##define a grid of the same size as the vector map
#r=raster(crs=proj4string(lower48),nrows=800,ncols=800)
#extent(r)=extent(lower48)
#r=as(r, 'SpatialGrid')
## indices not to plot, outside map
#idw.idx=!is.na(over(r,lower48)$iso_a2)
# read logo
#logo <- readPNG("~/Dropbox/radar/NEXRAD/fullcontinent/logo/CL_logo_stack_RGB_inv.png")
#res <- dim(logo)[1:2]
#save(lower48,lower48wgs,radarInfo,r,idw.idx,logo,res,file="basemap.RData")
# load basemap and radarInfo
load("basemap.RData")
# load variogram model for interpolation
load("vgModel.RData")

##########################################################
#  download files
##########################################################

setwd(VPDIR)
cat(paste("SCRIPT: making mosaic for",DATE,"\n"))
# select and load files
timer=system.time(keys<-do.call(rbind,lapply(radarInfo$radar,function(x) s3file(x,DATE))))
cat(paste("SCRIPT: selected", nrow(keys), "profiles for download in",round(timer["elapsed"],2),"seconds\n"))
# abort if too few radars available
stopifnot(nrow(keys)>RADARS_MIN)
timer=system.time(lapply(keys$Key,function(x) save_object(x,bucket=BUCKET,file=basename(x))))
cat(paste("SCRIPT: downloaded profiles to", VPDIR, "in",round(timer["elapsed"],2),"seconds\n"))
timer=system.time(vps<-suppressWarnings(readvp.list(list.files(pattern="*.h5"))))
cat(paste("SCRIPT: loaded",length(vps),"profiles in",round(timer["elapsed"],2),"seconds\n"))

##########################################################
#  interpolate and plot
##########################################################

# we interpolate with static variogram model, to avoid occasional singularities in variogram fit
vps.idw=ipol.vplist(vps,log=T,method="krige",variogram=vgModel)
# reset the date attribute to the requested value
attributes(vps.idw)$date=DATE
# plot to file
outputfile=paste(VPDIR,strftime(DATE,"/mosaic_%Y%m%d%H%M.jpg"),sep="")
plotidw(vps.idw,main="mtr",bg="black",linecol="white",zlim=c(50,50000),file=outputfile,closedev=F,legend.lab="Migration traffic rate [thousands/km/h]",variogram=vgModel)
add_barbs(vps)
# add day-night terminator
term=get_terminator(DATE)
if(!is.null(term)) points(term@coords[,1],term@coords[,2],col='yellow',type='l',lwd=2)
garbage=dev.off()

##########################################################
#  upload to S3
##########################################################

if(S3UPLOAD){
  cat(paste("uploading",outputfile),"...")
  success=put_object(outputfile,strftime(DATE,"mosaic/%Y/%m/%d/mosaic_%Y%m%d%H%M.jpg"),BUCKET,acl="public-read")
  # update the filenames text file
  filenames_local=paste(VPDIR,"/filenames.txt",sep="")
  # if after 22 UTC, we group the file with the next day
  if(as.numeric(strftime(DATE, format="%H"))>=HOUR_NEW_DAY){
    filenames_remote=strftime(DATE+24*3600,"mosaic/%Y/%m/%d/filenames.txt")
  } else filenames_remote=strftime(DATE,"mosaic/%Y/%m/%d/filenames.txt")
  # try to grab the filenames.txt from S3
  success=tryCatch(save_object(filenames_remote,BUCKET,filenames_local),error=function(x) FALSE)
  if(success==FALSE) cat(paste("starting new filenames.txt"))
  write(strftime(DATE,"%Y/%m/%d/mosaic_%Y%m%d%H%M.jpg"),file=filenames_local,append=TRUE)
  # put a copy filenames.txt down the date tree
  success=put_object(filenames_local,filenames_remote,BUCKET,acl="public-read")
  # also put a copy at the root
  success=put_object(filenames_local,"mosaic/filenames.txt",BUCKET,acl="public-read")
  cat("done\n")
}

#copy imagery to optional local directory
if(length(FIGDIR)!=0) file.copy(outputfile,paste(FIGDIR,"/",basename(outputfile),sep=""))

#clean up
unlink(VPDIR,recursive=TRUE)
