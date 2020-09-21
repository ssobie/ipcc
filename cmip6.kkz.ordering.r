##Order the set of CMIP5 models using climdex projex and KKZ algorithm

source('/storage/home/ssobie/code/repos/ipcc/scenario_selection/KKZ.R',chdir=TRUE)

standard.climdex.names <- c('cddETCCDI','csdiETCCDI','cwdETCCDI','dtrETCCDI',
                         'fdETCCDI','gslETCCDI','idETCCDI',
                         'prcptotETCCDI',
                         'r10mmETCCDI','r1mmETCCDI','r20mmETCCDI',
                         'r95pETCCDI','r99pETCCDI',
                         'rx1dayETCCDI','rx5dayETCCDI',
                         'sdiiETCCDI',
                         'suETCCDI',
                         'tn10pETCCDI','tn90pETCCDI',
                         'tnnETCCDI','tnxETCCDI',
                         'trETCCDI',
                         'tx10pETCCDI','tx90pETCCDI',
                         'txnETCCDI','txxETCCDI',
                         'wsdiETCCDI')



##-----------------------------------------------------------------------------------
region <- 'NWN'


read.dir <- '/storage/data/climate/CMIP6/KKZ/'
write.dir <- '/storage/data/climate/CMIP6/KKZ/model_rankings/'

method.files <- list.files(path=read.dir,pattern=paste0(region,'_'))
cmip6.files <- method.files[grep('2081-2100_1986-2005',method.files)]

projection.matrix <- matrix(NA,nrow=length(cmip6.files),ncol=27)

for (i in seq_along(cmip6.files)) {
   file <- cmip6.files[i]
   load(paste0(read.dir,file))
   ix <- names(clipped.anomaly.median) %in% standard.climdex.names  
   projection.matrix[i,] <- as.numeric(clipped.anomaly.median)[ix]

}

colnames(projection.matrix) <- standard.climdex.names
norm.proj.matrix <- round(scale(projection.matrix),5)
projection.matrix <- round(projection.matrix,5)

climdex.kkz <- subset.kkz(norm.proj.matrix, n.cases=nrow(norm.proj.matrix),
                                     newdata=projection.matrix)
models <- unlist(lapply(strsplit(cmip6.files[as.numeric(rownames(climdex.kkz$cases))],'_'),function(x){return(x[2])}))
runs <- unlist(lapply(strsplit(cmip6.files[as.numeric(rownames(climdex.kkz$cases))],'_'),function(x){return(x[4])}))

all.range.diff <- apply(apply(norm.proj.matrix, 2, range),2,diff)
all.range.thresh <- 0.9 * all.range.diff
model.thresh <- floor(0.9*length(models))


kkz.all <- climdex.kkz$cases.newdata[1,]
kkz.all <- rbind(kkz.all,climdex.kkz$cases.newdata[2,])

kkz.set <- climdex.kkz$cases[1,]
kkz.set <- rbind(kkz.set,climdex.kkz$cases[2,])

for(i in 3:nrow(climdex.kkz$cases.newdata)) {

    kkz.set <- rbind(kkz.set,climdex.kkz$cases[i,])
    kkz.all <- rbind(kkz.all,climdex.kkz$cases.newdata[i,])
    ##kkz.set.range <- apply(kkz.set, 2, range)
    kkz.set.mins <- ceiling(apply(kkz.set, 2, min)*1000)/1000
    kkz.set.maxs <- floor(apply(kkz.set, 2, max)*1000)/1000
    kkz.set.range <- rbind(kkz.set.mins,kkz.set.maxs)
    kkz.all.range <- apply(kkz.all, 2, range)
    kkz.set.diff <- apply(apply(kkz.set, 2, range),2,diff)
    print(paste0('Number of Models: ',i))
    print('Number of variables with 90% of range covered')
    print(paste0(sum(kkz.set.diff >= all.range.thresh),' of ',ncol(projection.matrix)))

    kkz.set.coverage <- colMeans(sweep(norm.proj.matrix, 2, kkz.set.range[1,], '>=') &
                                 sweep(norm.proj.matrix, 2, kkz.set.range[2,], '<='))
    kkz.lower.bound <- sum(apply(sweep(norm.proj.matrix, 2, kkz.set.range[1,], '>='),2,sum) >= model.thresh)
    print(paste0('KKZ lower bound: ',kkz.lower.bound))
    kkz.upper.bound <- sum(apply(sweep(norm.proj.matrix, 2, kkz.set.range[2,], '<='),2,sum) >= model.thresh)
    print(paste0('KKZ upper bound: ',kkz.upper.bound))

    kkz.all.coverage <- colMeans(sweep(projection.matrix, 2, kkz.all.range[1,], '>=') &
                                 sweep(projection.matrix, 2, kkz.all.range[2,], '<='))
    kkz.all.lower <- sum(apply(sweep(projection.matrix, 2, kkz.all.range[1,], '>='),2,sum) >= model.thresh)
    print(paste0('KKZ all lower bound: ',kkz.all.lower))
    kkz.all.upper <- sum(apply(sweep(projection.matrix, 2, kkz.all.range[2,], '<='),2,sum) >= model.thresh)
    print(paste0('KKZ all upper bound: ',kkz.all.upper))
    if (kkz.all.lower == 27 & kkz.all.upper == 27 ) { 
        print(cbind(models[1:i],runs[1:i]))
        output <- rbind(c('Order','Model','Run'),cbind(1:i,models[1:i],runs[1:i]))
        write.table(output,file=paste0(write.dir,region,'_ordered_CMIP6_climdex_KKZ.csv'),quote=F,row.name=F,col.name=F,sep=',')
        browser()
    }

    print('Model Added')

    ##print(paste0('KKZ Ordered Coverage ',kkz.set.coverage >= 0.9))
    print(paste0('KKZ Ordered Coverage ',mean(kkz.set.coverage >= 0.9)))
    if (sum(kkz.set.coverage >=0.9)==27) {
        print(models[1:i])
        browser()
    }

}