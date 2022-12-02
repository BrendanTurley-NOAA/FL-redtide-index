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
table(year(habs$date))

### modis derived parms
parm_out <- c('chlor_a','chl_anom','nflh','nflh_anom','rrs_667','abi','bbp_carder','bbp_morel','ssnlw488','rbd','kbbi','cm_bbp')

### load bathymetry and grid
setwd('~/Documents/nasa/data/lowres_4km')
load('bathy.RData')

lonlat_modis <- expand.grid(lon=lon_modis,lat=lat_modis)
lonlat_modis$lon_m <- cut(lonlat_modis$lon,vec_brk(lon_modis))
lonlat_modis$lat_m <- cut(lonlat_modis$lat,vec_brk(lat_modis))

### cut by bathymetry grid
habs$lon_m <- cut(habs$LONGITUDE,vec_brk(lon_modis))
habs$lat_m <- cut(habs$LATITUDE,vec_brk(lat_modis))

satsst_out <- data.frame(matrix(NA,nrow(habs),13))
names(satsst_out) <- c(parm_out,'sst')
habs_n <- cbind(habs,satsst_out)
yrs <- 2003:2021
pb <- txtProgressBar(min = 0, max = length(yrs), style = 3)
for(i in 1:length(yrs)){
  inx1 <- which(year(habs$date)==yrs[i])
  
  setwd('~/Documents/nasa/data/lowres_4km')
  modis <- readRDS(paste0('aqua_modis_daily_input_',yrs[i],'_4km.rds')) # upside down
  
  sst <- readRDS(paste0('mursst_daily_',yrs[i],'.rds')) # upside up
  sst <- sst[,dim(sst)[2]:1,] # flip tp make consistent
  
  for(j in inx1){
    tmp <- habs[j,]
    jday <- yday(habs$date[j])
    jday_1m <- jday-1
    jday_1p <- jday+1
    
    if(jday==1){
      jday_1m <- jday
      jday_1p <- jday+2
    }
    if(jday==365 | jday==366){
      jday_1m <- jday-2
      jday_1p <- jday
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
        sat_tmp <- modis[k,ind[1],ind[2],jday_1m:jday_1p]
        inx2 <- which(!is.na(sat_tmp))
        if(length(inx2>0)){
          sat_tmp <- modis[k,ind[1],ind[2],(jday_1m:jday_1p)[inx2[1]]]
        }
        if(all(is.na(sat_tmp))){
          ### find values within spaital window
          sat_tmp <- modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday]
          inx3 <- which(!is.na(sat_tmp[1,]))
          if(length(inx3>0)){
            sat_tmp <- mean(modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday],na.rm=T)
          }
          rm(inx3)
        }
        if(all(is.na(sat_tmp))){
          ### find values within spaital and temporal window
          sat_tmp <- modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday_1m:jday_1p]
          inx4 <- which(!is.na(sat_tmp[1,,]))
          if(length(inx4>0)){
            sat_tmp <- mean(modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday_1m:jday_1p],na.rm=T)
          }
          rm(inx4)
        }
        sat_return[k] <- sat_tmp
        rm(sat_tmp,inx2)
      }
    }
    
    habs_n[j,29:41] <- c(sat_return,sst_return)
    
    rm(tmp,jday,ind,ind_m,sat_return,sst_return)
  }
  rm(inx1,inx_na,modis,sst)
  setTxtProgressBar(pb, i)
}


habs_n$week <- week(habs_n$date)
names(habs_n)[c(3:4,6,10,26:42)]
habs_reduce <- habs_n[,c(3:4,6,10,26:42)]
habs_reduce$year <- year(habs_reduce$date)
habs_reduce$month <- month(habs_reduce$date)
habs_reduce$yday <- yday(habs_reduce$date)
habs_reduce$ygm <- paste(habs_reduce$year,habs_reduce$month,habs_reduce$week,habs_reduce$yday,habs_reduce$lon_m,habs_reduce$lat_m)
habs_agg <- aggregate(cbind(LATITUDE,LONGITUDE,SAMPLE_DEPTH,CELLCOUNT,date,chlor_a,chl_anom,nflh,nflh_anom,rrs_667,abi,bbp_carder,bbp_morel,ssnlw488,rbd,kbbi,cm_bbp,sst,year,month,yday,week)~ygm,data=habs_reduce,mean,na.rm=T)

habs_agg$date <- as.Date(habs_agg$yday,origin=paste0(habs_agg$year,'-01-01'))
habs_agg$pa100k <- rep(0,nrow(habs_agg))
habs_agg$pa100k[habs_agg$CELLCOUNT>=100000] <- 1


bathy <- data.frame(longitude=lonlat_modis$lon,
                    latitude=lonlat_modis$lat,
                    depth_m=as.vector(bathy_agg_modis[,ncol(bathy_agg_modis):1])) ### bathymetry is upside up; flipped to be consistent with other inputs
bathy$lon_m <- cut(bathy$longitude,vec_brk(lon_modis))
bathy$lat_m <- cut(bathy$latitude,vec_brk(lat_modis))

habs_agg$lon_m <- cut(habs_agg$LONGITUDE,vec_brk(lon_modis))
habs_agg$lat_m <- cut(habs_agg$LATITUDE,vec_brk(lat_modis))

hab_bathy <- merge(habs_agg,bathy[,-c(1,2)],by=c('lon_m','lat_m'),all.x=T) # don't include superfluous lon/lats
setwd('~/Documents/nasa/data/lowres_4km')
write.csv(hab_bathy,'habs_covariates.csv')

length(which(hab_bathy$depth_m>2))
hist(hab_bathy$depth_m[which(hab_bathy$depth_m>=0)])
plot(hab_bathy$LONGITUDE,hab_bathy$LATITUDE,cex=log(hab_bathy$depth_m),asp=1)
plot(hab_bathy$LONGITUDE,hab_bathy$LATITUDE,asp=1)

hab_bathy$depth_m[which(hab_bathy$depth_m>=0)] <- 0
plot(hab_bathy$LONGITUDE,hab_bathy$LATITUDE,cex=log(-hab_bathy$depth_m),asp=1)

test <- na.omit(habs_agg)
par(mfrow=c(2,1))
plot(habs_agg$date,habs_agg$CELLCOUNT+1,log='y')
plot(habs$date,habs$CELLCOUNT+1,log='y')

plot(habs_agg$date,habs_agg$pa100k)


table(habs_agg$year,habs_agg$month)
table(habs_agg$year,habs_agg$week)

habs$year <- year(habs$date)
habs$month <- month(habs$date)
habs$week <- week(habs$date)
table(habs$year,habs$month)
table(habs$year,habs$week)

names(habs_agg)

mod1 <- glm(pa100k~chlor_a,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~chl_anom,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~nflh,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~nflh_anom,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~rrs_667,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~abi,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~nflh_anom,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~bbp_carder,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~bbp_morel,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~ssnlw488,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~rbd,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~kbbi,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~cm_bbp,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~sst,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

mod1 <- glm(pa100k~yday,data=habs_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')


library(mgcv)

AllModel  <- gam(as.factor(pa100k) ~ as.factor(month) + te(chl_anom,week) + te(chlor_a,week) + te(bbp_carder,bbp_morel) + 
                   te(sst, depth_m) + te(rrs_667,week) + te(LONGITUDE,LATITUDE) + te(rbd,week) + te(depth_m), 
                 data=hab_bathy, family = binomial, select=TRUE, method="REML")
save(AllModel, file = "AllModel_initial.RData")
setwd('~/Documents/nasa/data/lowres_4km')
load('AllModel_initial.RData')
p <- plot(AllModel, pages=1, se=TRUE, cex.axis=2, cex.lab=1.5)
summary(AllModel)
