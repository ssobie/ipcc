##Script to plot time series of frost free days for Vancouver Intl.

library(ncdf4)
library(PCICt)
library(rgdal)
library(rgeos)
library(zoo)
library(scales)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)


##---------------------------------------------------------------------
##Read in GCM data

load_time_series <- function(var.name,region,scenario,read.dir) {

   files <- list.files(path=read.dir,pattern=paste0(var.name,'_',region))
   scen.file <- files[grep(scenario,files)] 

   time.series <- read.csv(paste0(read.dir,scen.file),header=T,as.is=T)

   return(time.series)

}

##---------------------------------------------------------------------
##Historical part of series

get_sub_series <- function(series,interval) {

   bounds <- strsplit(interval,'-')[[1]]

   dates <- as.Date(series[,1])
   yst <- head(grep(bounds[1],format(dates,'%Y')),1)
   yen <- tail(grep(bounds[2],format(dates,'%Y')),1)

   sub.series <- series[yst:yen,-1]
   sub.dates <- dates[yst:yen]
   rv <- list(series=sub.series,dates=sub.dates)
   return(rv)
}

##---------------------------------------------------------------------
##*********************************************************************

region <- 'canada_boundary'
cmip <- 'cmip6'
var.name <- 'tas'
period <- 'annual'

past.int <- switch(cmip,cmip5='1850-2005',cmip6='1850-2014')
proj.int <- switch(cmip,cmip5='2006-2100',cmip6='2015-2100')


base.dir <- '/storage/data/projects/rci/data/cas/'
read.dir <- paste0(base.dir,cmip,'/time_series/',region,'/')
plot.dir <- paste0(base.dir,cmip,'/plots/')

plot.title <- paste0('Annual Average Temperature Anomalies for Canada')
plot.title <- paste0('Temperature Change (Canada: Annual)')

plot.file <- paste0(plot.dir,region,'.average.',var.name,'.anomalies.',cmip,'.1971-2000.png')

##-----------------------------
##Observations
##GISTEMP
gis.file <- paste0("tas_gistemp_",region,"_time_series_anomalies_obs_1880-2017.csv")
gistemp <- read.csv(paste0(base.dir,"gridded_observations/time_series/",gis.file),header=TRUE,as.is=T)
gis.dates <- as.Date(gistemp[,1])

g.yst <- head(grep(1950,format(gis.dates,'%Y')),1)
g.yen <- tail(grep(2014,format(gis.dates,'%Y')),1)


##-----------------------------
tas.126 <- load_time_series(var.name,region,scenario='ssp126',read.dir)
tas.126.mean <- apply(tas.126[,-1],1,mean,na.rm=T)
tas.126.future <- get_sub_series(tas.126,proj.int)


tas.245 <- load_time_series(var.name,region,scenario='ssp245',read.dir)
tas.245.mean <- apply(tas.245[,-1],1,mean,na.rm=T)
tas.245.future <- get_sub_series(tas.245,proj.int)

tas.585 <- load_time_series(var.name,region,scenario='ssp585',read.dir)
tas.585.mean <- apply(tas.585[,-1],1,mean,na.rm=T)
tas.585.future <- get_sub_series(tas.585,proj.int)

tas.hist <- get_sub_series(tas.245,past.int)

x <- as.Date(tas.245[,1])

xtks <- pretty(x)
x.axis <- format(xtks,'%Y')

ylim <- c(-3.5,16)
y.axis <- pretty(ylim)


png(plot.file,width=6,height=3.5,units='in',res=300,pointsize=6,bg='white')
par(mar=c(5,5,5,3))
cx <- 1.2
plot(c(),xlim=c(as.Date('1850-01-01'),as.Date('2100-12-31')),
         ylim=ylim,yaxs='i',xaxs='i',
         main=plot.title,xlab='Date',ylab='Temperature (\u00B0C)',
         cex.axis=cx,cex.lab=cx,cex.main=cx+0.25,
         col.main='gray22',col.lab='gray22',axes=F)

axis(1,at=xtks,label=x.axis,cex.axis=cx,col='gray22',col.axis='gray22')
axis(2,at=y.axis,label=y.axis,cex.axis=cx,col='gray22',col.axis='gray22')

##Background
apply(tas.245[,-1],2,function(y,x){lines(x,y,col=alpha('gray',0.25))},x)
lines(x,apply(tas.245[,-1],1,mean,na.rm=T),lwd=3,col='gray15')

##Historical
apply(tas.hist$series,2,function(y,x){lines(x,y,col=alpha('gray',0.25))},tas.hist$dates)
lines(tas.hist$dates,apply(tas.hist$series,1,mean,na.rm=T),lwd=3,col='gray15')

##Future
apply(tas.585.future$series,2,function(y,x){lines(x,y,col=alpha('red',0.2))},tas.585.future$dates)
apply(tas.245.future$series,2,function(y,x){lines(x,y,col=alpha('orange',0.2))},tas.245.future$dates)
apply(tas.126.future$series,2,function(y,x){lines(x,y,col=alpha('blue',0.2))},tas.126.future$dates)

lines(tas.585.future$dates,apply(tas.585.future$series,1,mean,na.rm=T),lwd=3,col='red')
lines(tas.245.future$dates,apply(tas.245.future$series,1,mean,na.rm=T),lwd=3,col='orange')
lines(tas.245.future$dates,apply(tas.126.future$series,1,mean,na.rm=T),lwd=3,col='blue')

lines(gis.dates[g.yst:g.yen],gistemp[g.yst:g.yen,2],col='darkgreen',lwd=1)


#lines(1951:2100,rcp26.series,lwd=4,col='blue')
#lines(1951:2100,rcp45.series,lwd=4,col='orange')
#lines(1951:2100,rcp85.series,lwd=4,col='red')
abline(h=seq(-10,20,2.5),col='lightgray',lty=3,lwd=1)
legend('topleft',legend=c('SSP5 8.5','SSP2 4.5','SSP1 2.6','Historical','GISTEMP'),
        col=c('red','orange','blue','gray15','darkgreen'),
        box.col='gray22',text.col='gray22',cex=cx,pch=15,
        pt.cex=2.95, y.intersp=0.8, title.adj=0.2, xjust=0)
box(which='plot',col='gray22')
dev.off()



browser()
rcp26.series <- apply(rcp26.data,1,mean,na.rm=T)
rcp45.series <- apply(rcp45.data,1,mean,na.rm=T)
rcp85.series <- apply(rcp85.data,1,mean,na.rm=T)


hist.90 <- rollmean(apply(rcp85.data,1,quantile,0.9,na.rm=T),rx)[px]
hist.10 <- rollmean(apply(rcp85.data,1,quantile,0.1,na.rm=T),rx)[px]

rcp26.90 <- rollmean(apply(rcp26.data,1,quantile,0.9),rx)[ix]
rcp26.10 <- rollmean(apply(rcp26.data,1,quantile,0.1),rx)[ix]

rcp45.90 <- rollmean(apply(rcp45.data,1,quantile,0.9,na.rm=T),rx)[ix]
rcp45.10 <- rollmean(apply(rcp45.data,1,quantile,0.1,na.rm=T),rx)[ix]
rcp85.90 <- rollmean(apply(rcp85.data,1,quantile,0.9,na.rm=T),rx)[ix]
rcp85.10 <- rollmean(apply(rcp85.data,1,quantile,0.1,na.rm=T),rx)[ix]




##plot.file <- paste0(plot.dir,'nfld.annual.tas.smoothed.png')
##png(plot.file,width=900,height=900)
par(mar=c(5,5,5,3))
plot(1956:2095,rcp26.mean,type='l',lwd=4,col='green',ylim=c(0,10),
     main='Annual Average Temperatures',xlab='Year',ylab='TAS (degC)',
     cex.axis=2,cex.lab=2,cex.main=2.5,xlim=c(1950,2110))
apply(rcp85.data,2,function(y,x){lines(x,y,lwd=1,col=alpha('red',0.5))},1951:2100)
apply(rcp45.data,2,function(y,x){lines(x,y,lwd=1,col=alpha('orange',0.3))},1951:2100)
apply(rcp26.data,2,function(y,x){lines(x,y,lwd=1,col=alpha('blue',0.3))},1951:2100)

lines(1951:2100,rcp26.series,lwd=4,col='blue')
lines(1951:2100,rcp45.series,lwd=4,col='orange')
lines(1951:2100,rcp85.series,lwd=4,col='red')
abline(h=seq(0,20,5),col='gray',lty=3,lwd=3)
legend('topleft',legend=c('RCP8.5','RCP4.5','RCP2.6'),col=c('red','orange','green'),cex=2,pch=15)
box(which='plot')
##dev.off()

years <- 2081:2100
tas <- rcp26.series[131:150]
rcp.ext <- data.frame(years=years,tas=tas)
rcp.fit <- lm(tas~years,rcp.ext)
x <- 2081:2110
y.26 <- rcp.fit$coefficients[2]*x + rcp.fit$coefficients[1]
lines(x,y.26,col='blue')

years <- 2081:2100
tas <- rcp45.series[131:150]
rcp.ext <- data.frame(years=years,tas=tas)
rcp.fit <- lm(tas~years,rcp.ext)
x <- 2081:2110
y.45 <- rcp.fit$coefficients[2]*x + rcp.fit$coefficients[1]
lines(x,y.45,col='orange')

years <- 2081:2100
tas <- rcp85.series[131:150]
rcp.ext <- data.frame(years=years,tas=tas)
rcp.fit <- lm(tas~years,rcp.ext)
x <- 2081:2110
y.85 <- rcp.fit$coefficients[2]*x + rcp.fit$coefficients[1]
lines(x,y.85,col='red')

ts <- 1951:2100

st.10s <- grep(1986,ts)
en.10s <- grep(2016,ts)

st.20s <- grep(2011,ts)
en.20s <- grep(2040,ts)

st.50s <- grep(2041,ts)
en.50s <- grep(2070,ts)

st.80s <- grep(2071,ts)
en.80s <- grep(2100,ts)

col.names <- c('RCPs','2010s','2020s','2050s','2080s','2100')
row.names <- c('RCP2.6','RCP4.5','RCP8.5')

anoms.26 <- c(mean(rcp26.series[st.10s:en.10s]),
              mean(rcp26.series[st.20s:en.20s]),
              mean(rcp26.series[st.50s:en.50s]),
              mean(rcp26.series[st.80s:en.80s]),
              mean(y.26[11:30]))

anoms.45 <- c(mean(rcp45.series[st.10s:en.10s]),
              mean(rcp45.series[st.20s:en.20s]),
              mean(rcp45.series[st.50s:en.50s]),
              mean(rcp45.series[st.80s:en.80s]),
              mean(y.45[11:30]))
              
anoms.85 <- c(mean(rcp85.series[st.10s:en.10s]),
              mean(rcp85.series[st.20s:en.20s]),
              mean(rcp85.series[st.50s:en.50s]),
              mean(rcp85.series[st.80s:en.80s]),
              mean(y.85[11:30]))


bc.anoms <- cbind(row.names,round(rbind(anoms.26,anoms.45,anoms.85),1))
colnames(bc.anoms) <- col.names

write.table(bc.anoms,file=paste0(plot.dir,'BC.gcm.bccaq2.tas.anomalies.csv'),quote=F,sep=',',col.names=T,row.names=F)
