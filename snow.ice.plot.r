##Script to make a multiplot figure of Snow and Ice data in the Canadian Arctic
library(scales)
ice.plot <- function(tas,ice.data) {

newx <- seq(-100, 400, by=0.5)
ice.fit <- matrix(nrow=length(newx),ncol=glen)

  tas.anoms <- vector(mode='list',length=glen)
  ice.anoms <- vector(mode='list',length=glen)

  for (i in 1:glen) {
     tas.mean <- mean(tas[c(11,12,13),i])
     tas.anom <- tas[,i]-tas.mean         
     tas.anoms[[i]] <- tas.anom
     ix <- order(tas.anom)

     ice.mean <- mean(ice.data[c(11,12,13),i])
     ice.anom <- (ice.data[,i] - ice.mean) ##/ice.mean
     ice.anoms[[i]] <- ice.anom
     ##points(ice.data[ix,i],tas.anom[ix],pch=16,col='gray')
     ##points(ice.anom[ix],tas.anom[ix],pch=16,col='gray')
     ##lines(ice.data[ix,i],tas.anom[ix])

     ##x <- as.vector(ice.data[ix,i])
     x <- as.vector(ice.anom[ix])
     y <- as.vector(tas.anom[ix])
     df <- data.frame(x = x,
                 y = y)
     mod <- lm(y ~ x, data = df)
     ##abline(mod,col='red')
     preds <- predict(mod, newdata = data.frame(x=newx))                     
     ice.fit[,i] <- preds
  }

  for (i in 1:20) {
      print(gcms[i])
      points(ice.anoms[[i]],tas.anoms[[i]],pch=i)
  }

  browser()
  test <- chull(unlist(ice.anoms),unlist(tas.anoms))
  hpts <- c(test,test[1])
  polygon(unlist(ice.anoms)[hpts],unlist(tas.anoms)[hpts],col=alpha('red',0.3),border=alpha('red',0.3))

  ice.lower <- apply(ice.fit,1,quantile,0.05)
  ice.upper <- apply(ice.fit,1,quantile,0.95)
  ice.mean <- apply(ice.fit,1,mean)

  lines(newx,ice.lower,col='blue',lwd=3)
  lines(newx,ice.mean,col='red',lwd=3)
  lines(newx,ice.upper,col='blue',lwd=3)
}


proj.dir <- '/storage/data/projects/rci/data/nrcan/'

read.dir <- paste0(proj.dir,'ice_snow_data/')
plot.dir <- paste0(proj.dir,'plots/')

arctic.temps <- read.csv(paste0(read.dir,'arctic_temps.csv'),as.is=T,header=F)
ice.free.nep <- read.csv(paste0(read.dir,'ice_free_nep.csv'),as.is=T,header=F)
ice.free.caa <- read.csv(paste0(read.dir,'ice_free_caa.csv'),as.is=T,header=F)

snow.cover.ea <- read.csv(paste0(read.dir,'snow_cover_eurasia.csv'),as.is=T,header=F)
snow.cover.na <- read.csv(paste0(read.dir,'snow_cover_north_america.csv'),as.is=T,header=F)


years <- arctic.temps[3:27,1]
models <- arctic.temps[1,-1]
runs <- arctic.temps[2,-1]

gcms <- paste0(models,'-r',runs)
csiro.mask <- grep('CSIRO',gcms)+1
glen <- length(gcms) - length(csiro.mask)
mask <- -1*c(1,csiro.mask)

tas.data <- arctic.temps[3:27,mask]
ice.nep.data <- ice.free.nep[3:27,mask]
ice.caa.data <- ice.free.caa[3:27,mask]

snow.ea.data <- snow.cover.ea[3:27,mask]
snow.na.data <- snow.cover.na[3:27,mask]


tas <- round(apply(tas.data,2,as.numeric),1)-273
ice.nep <- round(apply(ice.nep.data,2,as.numeric),1)
ice.caa <- round(apply(ice.caa.data,2,as.numeric),1)

snow.ea <- round(apply(snow.ea.data,2,as.numeric),1)
snow.na <- round(apply(snow.na.data,2,as.numeric),1)

if (1==0) {
  plot.file <- paste0(plot.dir,'arctic_ice_free_days_and_snow_cover3.png')
  png(plot.file,width=1200,height=1200)
  par(mar=c(5,5,4,2))
  par(mfrow=c(2,2))

  plot(c(),xlim=c(-50,300),ylim=c(-2,16),cex=2,main='Northeast Passage',
                   xlab='Ice Free Day Anomaly',ylab='Arctic Temperature Anomalies (\u00B0C)',
                   cex.lab=2,cex.axis=2,cex.main=2)
  ice.plot(tas,ice.nep)
  legend('topleft',leg=c('RCP2.6','RCP4.5','RCP8.5'),
          col=c(alpha('blue',0.3),alpha('green',0.3),alpha('red',0.3)),
          pch=c(15,15,15),cex=2)

  ## text(x=50,y=11,'North East\n Passage',cex=1.5) 
  ##text(x=-1,y=11,'North East\n Passage',cex=2)
  box(which='plot')

  plot(c(),xlim=c(-50,300),ylim=c(-2,16),cex=2,main='Canadian Arctic Archipelago',
                   xlab='Ice Free Day Anomaly',ylab='Arctic Temperature Anomalies (\u00B0C)',
                   cex.lab=2,cex.axis=2,cex.main=2)
  ice.plot(tas,ice.caa)
  ##text(x=65,y=11,'Canadian Arctic\nArchipelago',cex=1.5)
  ##text(x=-1,y=11,'Canadian Arctic\nArchipelago',cex=2)
  box(which='plot')
}
  plot(c(),xlim=c(-90,25),ylim=c(-2,16),cex=2,main='North America',
                   xlab='Snow Cover Day Anomaly',ylab='Arctic Temperature Anomalies (\u00B0C)',
                   cex.lab=2,cex.axis=2,cex.main=2)
  ice.plot(tas,snow.na)
  ##text(x=300,y=11,'North America',cex=1.5)
  ##text(x=0.12,y=11,'North America',cex=2)
  box(which='plot')

if (1==0) {
  plot(c(),xlim=c(-90,25),ylim=c(-2,16),cex=2,main='Eurasia',
                   xlab='Snow Cover Day Anomaly',ylab='Arctic Temperature Anomalies (\u00B0C)',
                   cex.lab=2,cex.axis=2,cex.main=2)
  ice.plot(tas,snow.ea)
  ##text(x=300,y=11,'Eurasia',cex=1.5)
  ##text(x=0.16,y=11,'Eurasia',cex=2)
  box(which='plot')


  dev.off()
}


# predicts + interval

# plot
##plot(y ~ x, data = df, type = 'n')
# add fill
##polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)
# model
##abline(mod)
# intervals
##lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
##lines(newx, preds[ ,2], lty = 'dashed', col = 'red')