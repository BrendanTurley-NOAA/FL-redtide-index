gc()

library(fields)
library(lubridate)

vec_brk <- function(input) {
  core <- (diff(input)/2) + input[1:(length(input)-1)]
  beg <- input[1] - (diff(input[1:2])/2)
  last <- input[length(input)] + (diff(input[(length(input)-1):(length(input))])/2)
  output <- c(beg,core,last)
  output
}

ind2sub <- function(ind,a){
  m <- nrow(a)
  r <- ((ind-1) %% m) + 1
  c <- floor((ind-1) / m) + 1
  return(c(r,c))
}

### load habs
setwd('~/Desktop/professional/projects/Postdoc_FL/data/habs')
habs <- read.csv('habsos_subset_FL03-21.csv')
habs$date <- ymd_hms(habs$date)
habs <- habs[order(habs$date),]
# habs$year <- year(habs$date)
# habs$month <- month(habs$date)
# habs$week <- week(habs$date)
# table(habs$year)
# table(habs$year,habs$month)
# table(habs$year,habs$week)
### proportion positive
length(which(habs$CELLCOUNT>1e5))/nrow(habs)

### modis derived parms
parm_out <- c('chlor_a','chl_anom','nflh','nflh_anom','rrs_667','ssnlw488','carder_bbp','morel_bbp','cm_bbp','abi','rbd','kbbi')

### load bathymetry and grid
setwd('~/Documents/nasa/data/lowres_4km')
load('bathy.RData')

lonlat_modis <- expand.grid(lon=lon_modis,lat=lat_modis)
lonlat_modis$lon_m <- cut(lonlat_modis$lon,vec_brk(lon_modis))
lonlat_modis$lat_m <- cut(lonlat_modis$lat,vec_brk(lat_modis))

bathy <- data.frame(longitude=lonlat_modis$lon,
                    latitude=lonlat_modis$lat,
                    depth_m=as.vector(bathy_agg_modis[,ncol(bathy_agg_modis):1])) ### bathymetry is upside up; flipped to be consistent with other inputs
bathy$lon_m <- cut(bathy$longitude,vec_brk(lon_modis))
bathy$lat_m <- cut(bathy$latitude,vec_brk(lat_modis))

### cut by bathymetry grid
habs$lon_m <- cut(habs$LONGITUDE,vec_brk(lon_modis))
habs$lat_m <- cut(habs$LATITUDE,vec_brk(lat_modis))

satsst_out <- data.frame(matrix(NA,nrow(habs),13))
names(satsst_out) <- c(parm_out,'sst')
habs_n <- cbind(habs,satsst_out)
yrs <- 2003:2021
pb <- txtProgressBar(min = 0, max = length(yrs), style = 3)
t1 <- system.time(
  for(i in 1:length(yrs)){
    inx1 <- which(year(habs$date)==yrs[i])
    
    # setwd('~/Documents/nasa/data/lowres_4km')
    setwd("/Volumes/Clovis/data/nasa/lowres_4km")
    modis <- readRDS(paste0('aqua_modis_daily_input_',yrs[i],'_4km.rds')) # upside down
    # if(yrs[i]>2003){ #
    #   modis_m <- readRDS(paste0('aqua_modis_daily_input_',yrs[i]-1,'_4km.rds')) # upside down
    #   modis_m <- modis_m[,,,(dim(modis_p)[4]-4):(dim(modis_p)[4])]
    # }
    # if(yrs[i]<2021){
    #   modis_p <- readRDS(paste0('aqua_modis_daily_input_',yrs[i]+1,'_4km.rds')) # upside down
    #   modis_p <- modis_p[,,,1:5]  
    # }
    
    sst <- readRDS(paste0('mursst_daily_',yrs[i],'.rds')) # upside up
    sst <- sst[,dim(sst)[2]:1,] # flip tp make consistent
    
    for(j in inx1){
      tmp <- habs[j,]
      jday <- yday(habs$date[j])
      jday_1m <- jday - 2 # (2023/03/28) increase window to 5 from 3
      jday_1p <- jday + 2 # (2023/03/28) increase window to 5 from 3
      
      if(jday_1m<1){ # would be more accurate to use preceding/following years data
        jday_1m <- 1
        jday_1p <- 5
      }
      if(jday_1p>365){
        jday_1m <- 360
        jday_1p <- 365
      }
      if(jday_1p>365 & leap_year(yrs[i])){
        jday_1m <- 361
        jday_1p <- 366
      }
      
      ind_m <- which(tmp$lon_m==lonlat_modis$lon_m & tmp$lat_m==lonlat_modis$lat_m)
      ind <- ind2sub(ind_m,matrix(NA,158,158))
      
      ### find values at exact time / location
      sat_return <- modis[,ind[1],ind[2],jday]
      sst_return <- sst[ind[1],ind[2],jday]
      
      inx_na <- which(is.na(sat_return))
      if(length(inx_na)>0){
        for(k in inx_na){
          ### find values within temporal window
          jwin <- jday_1m:jday_1p
          sat_tmp <- modis[k,ind[1],ind[2],jwin]
          # if(yr[i]>2003 & jday<3){ # would be more accurate to use preceding/following years data
          #   sat_tmp <- c(modis_m[k,ind[1],ind[2],(5+(jwin)[which((jwin)<1)])],
          #                modis[k,ind[1],ind[2],(jwin)[which((jwin)>0)]])
          #   add <- length(which((jwin)<1))
          #   jwin <- jwin+add
          #   jday <- jday+add
          # }
          # if(yr[i]<2021 & jday>364){
          #   sat_tmp <- c(modis[k,ind[1],ind[2],(jwin)[which((jwin)<366)]],
          #                modis_m[k,ind[1],ind[2],((jwin)[which((jwin)>365)]-365)])
          #   add <- length(which((jwin)>366))
          #   jwin <- jwin-add
          #   jday <- jday-add
          # }
          inx2 <- which(!is.na(sat_tmp))
          if(length(inx2)>0){
            # sat_tmp <- modis[k,ind[1],ind[2],(jwin)[inx2[1]]] # takes first index not NA
            nearest_sat <- inx2[which.min(abs((jwin)[inx2]-jday))]
            sat_tmp <- sat_tmp[nearest_sat] # takes temporally nearest satellite data
          }
          #### --- (2023/03/28) this chunk commented out for 4km model to reduce noise and increase accuracy of low res model
          # if(all(is.na(sat_tmp))){
          #   ### find values within spaital window
          #   sat_tmp <- modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday]
          #   inx3 <- which(!is.na(sat_tmp[1,]))
          #   if(length(inx3)>0){
          #     sat_tmp <- mean(modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday],na.rm=T)
          #   }
          #   rm(inx3)
          # }
          # if(all(is.na(sat_tmp))){
          #   ### find values within spaital and temporal window
          #   sat_tmp <- modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday_1m:jday_1p]
          #   inx4 <- which(!is.na(sat_tmp[1,,]))
          #   if(length(inx4)>0){
          #     sat_tmp <- mean(modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday_1m:jday_1p],na.rm=T)
          #   }
          #   rm(inx4)
          # }
          #### --- (2023/03/28) this chunk commented out 
          sat_return[k] <- ifelse(length(sat_tmp)==1, sat_tmp, mean(sat_tmp, na.rm=T))
          rm(sat_tmp,inx2)
        }
      }
      
      if(is.na(sst_return)){
        ### find values within temporal window
        jwin <- jday_1m:jday_1p
        sst_tmp <- sst[ind[1],ind[2],jwin]
        inx3 <- which(!is.na(sst_tmp))
        if(length(inx3)>0){
          # sat_tmp <- modis[k,ind[1],ind[2],(jday_1m:jday_1p)[inx3[1]]] # takes first index not NA
          nearest_sat <- inx3[which.min(abs((jwin)[inx3]-jday))]
          sst_tmp <- sst_tmp[nearest_sat] # takes temporally nearest satellite data
        }
        sst_return <- ifelse(length(sst_tmp)==1, sst_tmp, mean(sst_tmp, na.rm=T))
        rm(sst_tmp,inx3)
      }
  
      
      habs_n[j,29:41] <- c(sat_return,sst_return)
      
      rm(tmp,jday,ind,ind_m,sat_return,sst_return)
    }
    rm(inx1,inx_na,modis,sst)
    setTxtProgressBar(pb, i)
  }
)
t1
# user   system  elapsed 
# 624.090  444.737 1205.007 
# user   system  elapsed 
# 674.875  453.739 1338.307
# user  system elapsed 
# 522.949 364.737 941.802 
setwd('~/Documents/nasa/data/lowres_4km')
# write.csv(habs_n,'habs_covariates_full_v2.csv',row.names = F) # (2023/03/23) v2 no longer uses a spatial window if satellite values NA and increases temporal window from 3 to 5; also adds SST window
habs_n <- read.csv('habs_covariates_full_v2.csv')
habs_n$date <- ymd_hms(habs_n$date)

names(habs_n)[c(3:4,6,10,26:41)]
habs_reduce <- habs_n[,c(3:4,6,10,26:41)]
habs_reduce$year <- year(habs_reduce$date)
habs_reduce$month <- month(habs_reduce$date)
habs_reduce$week <- week(habs_reduce$date)
habs_reduce$yday <- yday(habs_reduce$date)
### weekly aggregates
habs_reduce$ygm <- paste(habs_reduce$year,habs_reduce$month,habs_reduce$week,habs_reduce$lon_m,habs_reduce$lat_m)
habs_agg <- aggregate(cbind(LATITUDE,LONGITUDE,SAMPLE_DEPTH,CELLCOUNT,date,chlor_a,chl_anom,nflh,nflh_anom,rrs_667,ssnlw488,carder_bbp,morel_bbp,cm_bbp,abi,rbd,kbbi,sst,year,month,yday,week)~ygm,
                      data=habs_reduce,mean,na.rm=T)
habs_agg$date <- as.Date(habs_agg$yday-1,origin=paste0(habs_agg$year,'-01-01'))
### create binary classifier
habs_agg$pa100k <- ifelse(habs_agg$CELLCOUNT>=1e5,1,0)
### assign to grid
habs_agg$lon_m <- cut(habs_agg$LONGITUDE,vec_brk(lon_modis))
habs_agg$lat_m <- cut(habs_agg$LATITUDE,vec_brk(lat_modis))
### add bathymetry
habs_covar_agg <- merge(habs_agg,bathy[,-c(1,2)],by=c('lon_m','lat_m'),all.x=T) # don't include superfluous lon/lats
### remove shallow water samples
# habs_covar_agg <- habs_covar_agg[-which(habs_covar_agg$depth_m>2),] # remove samples taken at altitude greater than 2 m
# updated 2023/03/28 to remove shallow water samples; remove samples taken at locations with depth greater than 200 m
habs_covar_agg1 <- habs_covar_agg[-which(habs_covar_agg$depth_m>(-1) |
                                           habs_covar_agg$depth_m<(-200)),]
setwd('~/Documents/nasa/data/lowres_4km')
write.csv(habs_covar_agg1,'habs_covariates_agg_v2.csv',row.names = F)

habs_covar_agg2 <- habs_covar_agg[-which(habs_covar_agg$depth_m>(-2) |
                                           habs_covar_agg$depth_m<(-200)),]
setwd('~/Documents/nasa/data/lowres_4km')
write.csv(habs_covar_agg2,'habs_covariates_agg_v2_2m.csv',row.names = F)

habs_covar_agg3 <- habs_covar_agg[-which(habs_covar_agg$depth_m>(-3) |
                                           habs_covar_agg$depth_m<(-200)),]
setwd('~/Documents/nasa/data/lowres_4km')
write.csv(habs_covar_agg3,'habs_covariates_agg_v2_3m.csv',row.names = F)

habs_covar_agg4 <- habs_covar_agg[-which(habs_covar_agg$depth_m>(-5) |
                                           habs_covar_agg$depth_m<(-200)),]
setwd('~/Documents/nasa/data/lowres_4km')
write.csv(habs_covar_agg4,'habs_covariates_agg_v2_5m.csv',row.names = F)

plot(habs_covar_agg1$LONGITUDE,habs_covar_agg1$LATITUDE,asp=1,pch='.',col=2)
points(habs_covar_agg2$LONGITUDE,habs_covar_agg2$LATITUDE,asp=1,pch='.',col=3)
points(habs_covar_agg3$LONGITUDE,habs_covar_agg3$LATITUDE,asp=1,pch='.',col=4)
points(habs_covar_agg4$LONGITUDE,habs_covar_agg4$LATITUDE,asp=1,pch='.')



test <- na.omit(habs_covar_agg)
dim(test)==dim(habs_covar_agg) ### should be true

par(mfrow=c(2,1))
plot(habs_covar_agg$date,habs_covar_agg$CELLCOUNT+1,log='y')
plot(habs$date,habs$CELLCOUNT+1,log='y')


### daily aggregates
habs_reduce$ygd <- paste(habs_reduce$year,habs_reduce$yday,habs_reduce$lon_m,habs_reduce$lat_m)
habs_agg2 <- aggregate(cbind(LATITUDE,LONGITUDE,SAMPLE_DEPTH,CELLCOUNT,date,chlor_a,chl_anom,nflh,nflh_anom,rrs_667,ssnlw488,carder_bbp,morel_bbp,cm_bbp,abi,rbd,kbbi,sst,year,month,yday,week)~ygd,
                      data=habs_reduce,mean,na.rm=T)
habs_agg2$date <- as.Date(habs_agg2$yday-1,origin=paste0(habs_agg2$year,'-01-01'))
### create binary classifier
habs_agg2$pa100k <- ifelse(habs_agg2$CELLCOUNT>=1e5,1,0)
### assign to grid
habs_agg2$lon_m <- cut(habs_agg2$LONGITUDE,vec_brk(lon_modis))
habs_agg2$lat_m <- cut(habs_agg2$LATITUDE,vec_brk(lat_modis))
### add bathymetry
habs_covar_agg1.2 <- merge(habs_agg2,bathy[,-c(1,2)],by=c('lon_m','lat_m'),all.x=T) # don't include superfluous lon/lats
### remove shallow water samples
# habs_covar_agg <- habs_covar_agg[-which(habs_covar_agg$depth_m>2),] # remove samples taken at altitude greater than 2 m
# updated 2023/03/28 to remove shallow water samples; remove samples taken at locations with depth greater than 200 m
habs_covar_agg1.3 <- habs_covar_agg1.2[-which(habs_covar_agg1.2$depth_m>(-1) |
                                             habs_covar_agg1.2$depth_m<(-200)),]
habs_covar_agg2.2 <- habs_covar_agg1.2[-which(habs_covar_agg1.2$depth_m>(-2) |
                                             habs_covar_agg1.2$depth_m<(-200)),]
habs_covar_agg3.2 <- habs_covar_agg1.2[-which(habs_covar_agg1.2$depth_m>(-3) |
                                             habs_covar_agg1.2$depth_m<(-200)),]
habs_covar_agg4.2 <- habs_covar_agg1.2[-which(habs_covar_agg1.2$depth_m>(-5) |
                                             habs_covar_agg1.2$depth_m<(-200)),]

### alternative; does change pa100k values
habs_agg2.1 <- aggregate(cbind(LATITUDE,LONGITUDE,SAMPLE_DEPTH,date,chlor_a,chl_anom,nflh,nflh_anom,rrs_667,ssnlw488,carder_bbp,morel_bbp,cm_bbp,abi,rbd,kbbi,sst,year,month,yday,week)~ygm,
                         data=habs_reduce,mean,na.rm=T)
habs_agg2.2 <- aggregate(CELLCOUNT~ygm,data=habs_reduce,max,na.rm=T)
habs_agg2.2 <- aggregate(CELLCOUNT~ygm,data=habs_reduce,quantile,.9,na.rm=T)
habs_agg2.2 <- aggregate(CELLCOUNT~ygm,data=habs_reduce,
                         function(x){if(length(x)<3){mean(x,na.rm=T)}else{mean(x[which(x>quantile(x,.9,na.rm=T))],na.rm=T)}})
habs_agg2.3 <- aggregate(CELLCOUNT~ygm,data=habs_reduce,length)
habs_agg2 <- merge(habs_agg2.1,habs_agg2.2,by=c('ygm'),all.x=T)
habs_agg2$date <- as.Date(habs_agg2$yday-1,origin=paste0(habs_agg2$year,'-01-01'))
# habs_aggt <- merge(habs_agg,habs_agg2,by=c('ygm'),all=T)
### create binary classifier
habs_agg2$pa100k <- ifelse(habs_agg2$CELLCOUNT>=1e5,1,0)
### assign to grid
habs_agg2$lon_m <- cut(habs_agg2$LONGITUDE,vec_brk(lon_modis))
habs_agg2$lat_m <- cut(habs_agg2$LATITUDE,vec_brk(lat_modis))
### add bathymetry
habs_covar_agg2 <- merge(habs_agg2,bathy[,-c(1,2)],by=c('lon_m','lat_m'),all.x=T) # don't include superfluous lon/lats
habs_covar_agg2 <- habs_covar_agg2[-which(habs_covar_agg2$depth_m>2),] # remove samples taken at altitude greater than 2 m
habs_covar_agg2 <- habs_covar_agg2[-which(habs_covar_agg2$depth_m<(-200)),] # remove samples taken at locations with depth greater than 200 m

setwd('~/Documents/nasa/data/lowres_4km')
# write.csv(habs_covar_agg2,'habs_covariates_agg2.csv',row.names = F)

test <- na.omit(habs_covar_agg2)
dim(test)==dim(habs_covar_agg2) ### should be true

par(mfrow=c(3,1))
plot(habs_covar_agg$date,habs_covar_agg$CELLCOUNT+1,log='y')
plot(habs_covar_agg2$date,habs_covar_agg2$CELLCOUNT+1,log='y')
plot(habs$date,habs$CELLCOUNT+1,log='y')

length(which(habs_covar_agg$pa100k==1))
length(which(habs_covar_agg2$pa100k==1))

