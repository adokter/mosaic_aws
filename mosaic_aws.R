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
setwd("/opt")

Sys.setenv(TZ="UTC")

VPDIR=paste(tempdir(),"/vp",sep="")
suppressWarnings(dir.create(VPDIR))

if(file.exists("~/.aws/credentials")){
  getkeys=function() read.csv("~/.aws/credentials",sep="=",stringsAsFactors = F,col.names = "key",strip.white=T)
  Sys.setenv(AWS_ACCESS_KEY_ID=getkeys()["aws_access_key_id",])
  Sys.setenv(AWS_SECRET_ACCESS_KEY=getkeys()["aws_secret_access_key",])
}

#Sys.setenv(MOSAIC_DATE=as.character("2017-10-20 01:00"))
DATE=as.POSIXct(if(Sys.getenv("MOSAIC_DATE")=="") as.character(Sys.time()) else Sys.getenv("MOSAIC_DATE"))

evening.palette <- colorRampPalette(c("#0a0a0a","#1775cb","white"),bias=1.0)(100)

OUTPUTDIR="~/Dropbox/radar/aws_visualisations/realtime"

# minimum number of radars for outputting an image
RADARS_MIN=120

##########################################################
#  functions
##########################################################

ipol.vplist = function(vplist,alpha=NA,log=F,alt.min=0, alt.max=Inf,quantity="mtr",method="idw"){
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
    vario.fit <- fit.variogram(vario.mtr, model=vgm(model='Sph'))
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
  text(usr[1] + xwidth/23, usr[3] + yheight/5, my.date.lab, col=linecol,cex = 5.5, font=2,pos=4)

  # plot logo
  if(!is.na(file)){
    rasterImage(logo, usr[2]+usr[2]/5.5, usr[3] + (yheight * 0.1) + yheight/20,
                usr[2] +usr[2]/5.5 - (((usr[3]) - (usr[3] - (yheight * 0.1)))*(res[2]/res[1])),
                usr[3] + yheight/20, interpolate=FALSE,angle=180)
  }

  bins=seq(zlim[1],zlim[2],length.out=length(col)+1)

  # plot radar locations
  plot(attributes(idw)$data.vp,add=T,col=linecol,lwd=0.1,pch='.')

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
                            cex.axis=5,pos=1,mgp = c(0, 3, 0),
                            lwd.ticks=2),
             legend.args=list(text=legend.lab, side=1, cex=.5, col=linecol, line=8)
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

s3file=function(radar,date,bucket="vol2bird",prefix="output",dt=900){
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
load("basemap.RData")

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
timer=system.time(lapply(keys$Key,function(x) save_object(x,bucket="vol2bird",file=basename(x))))
cat(paste("SCRIPT: downloaded profiles to", VPDIR, "in",round(timer["elapsed"],2),"seconds\n"))
timer=system.time(vps<-suppressWarnings(readvp.list(list.files(pattern="*.h5"))))
cat(paste("SCRIPT: loaded",length(vps),"profiles in",round(timer["elapsed"],2),"seconds\n"))

##########################################################
#  interpolate
##########################################################
vps.idw=ipol.vplist(vps,log=T,method="krige")
# reset the date attribute to the requested value
attributes(vps.idw)$date=DATE
# plot to file
outputfile=paste(VPDIR,strftime(DATE,"/mosaic_%Y%m%d%H%M.jpg"),sep="")
plotidw(vps.idw,main="mtr",bg="black",linecol="white",zlim=c(50,50000),file=outputfile,closedev=F,legend.lab="Migration traffic rate [kbird/km/h]")
add_barbs(vps)
garbage=dev.off()
# upload to S3
cat(paste("uploading",outputfile),"...")
success=put_object(outputfile,strftime(DATE,"mosaic/%Y/%m/%d/mosaic_%Y%m%d%H%M.jpg"),"vol2bird",acl="public-read")
cat("done\n")