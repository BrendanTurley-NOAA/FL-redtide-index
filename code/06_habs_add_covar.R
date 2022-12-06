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

### modis derived parms
parm_out <- c('chlor_a','chl_anom','nflh','nflh_anom','rrs_667','abi','bbp_carder','bbp_morel','ssnlw488','rbd','kbbi','cm_bbp')

### load bathymetry and grid
setwd('~/Documents/nasa/data/lowres_4km')
load('bathy.RData')

bathy <- data.frame(longitude=lonlat_modis$lon,
                    latitude=lonlat_modis$lat,
                    depth_m=as.vector(bathy_agg_modis[,ncol(bathy_agg_modis):1])) ### bathymetry is upside up; flipped to be consistent with other inputs
bathy$lon_m <- cut(bathy$longitude,vec_brk(lon_modis))
bathy$lat_m <- cut(bathy$latitude,vec_brk(lat_modis))

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
t1 <- system.time(
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
          if(length(inx2)>0){
            sat_tmp <- modis[k,ind[1],ind[2],(jday_1m:jday_1p)[inx2[1]]]
          }
          if(all(is.na(sat_tmp))){
            ### find values within spaital window
            sat_tmp <- modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday]
            inx3 <- which(!is.na(sat_tmp[1,]))
            if(length(inx3)>0){
              sat_tmp <- mean(modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday],na.rm=T)
            }
            rm(inx3)
          }
          if(all(is.na(sat_tmp))){
            ### find values within spaital and temporal window
            sat_tmp <- modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday_1m:jday_1p]
            inx4 <- which(!is.na(sat_tmp[1,,]))
            if(length(inx4)>0){
              sat_tmp <- mean(modis[k,(ind[1]-1):(ind[1]+1),(ind[2]-1):(ind[2]+1),jday_1m:jday_1p],na.rm=T)
            }
            rm(inx4)
          }
          sat_return[k] <- ifelse(length(sat_tmp)==1,sat_tmp,NA)
          rm(sat_tmp,inx2)
        }
      }
      
      habs_n[j,29:41] <- c(sat_return,sst_return)
      
      rm(tmp,jday,ind,ind_m,sat_return,sst_return)
    }
    rm(inx1,inx_na,modis,sst)
    setTxtProgressBar(pb, i)
  }
)
# user   system  elapsed 
# 674.875  453.739 1338.307
# write.csv(habs_n,'habs_covariates_full.csv',row.names = F)

names(habs_n)[c(3:4,6,10,26:41)]
habs_reduce <- habs_n[,c(3:4,6,10,26:41)]
habs_reduce$year <- year(habs_reduce$date)
habs_reduce$month <- month(habs_reduce$date)
habs_reduce$week <- week(habs_reduce$date)
habs_reduce$yday <- yday(habs_reduce$date)
### weekly aggregates
habs_reduce$ygm <- paste(habs_reduce$year,habs_reduce$month,habs_reduce$week,habs_reduce$lon_m,habs_reduce$lat_m)
habs_agg <- aggregate(cbind(LATITUDE,LONGITUDE,SAMPLE_DEPTH,CELLCOUNT,date,chlor_a,chl_anom,nflh,nflh_anom,rrs_667,abi,bbp_carder,bbp_morel,ssnlw488,rbd,kbbi,cm_bbp,sst,year,month,yday,week)~ygm,data=habs_reduce,mean,na.rm=T)
habs_agg$date <- as.Date(habs_agg$yday-1,origin=paste0(habs_agg$year,'-01-01'))
### create binary classifier
habs_agg$pa100k <- rep(0,nrow(habs_agg))
habs_agg$pa100k[habs_agg$CELLCOUNT>=100000] <- 1
### assign to grid
habs_agg$lon_m <- cut(habs_agg$LONGITUDE,vec_brk(lon_modis))
habs_agg$lat_m <- cut(habs_agg$LATITUDE,vec_brk(lat_modis))
### add bathymetry
habs_covar_agg <- merge(habs_agg,bathy[,-c(1,2)],by=c('lon_m','lat_m'),all.x=T) # don't include superfluous lon/lats
habs_covar_agg <- habs_covar_agg[-which(habs_covar_agg$depth_m>2),] # remove samples taken at altitude greater than 2 m
habs_covar_agg <- habs_covar_agg[-which(habs_covar_agg$depth_m<(-200)),] # remove samples taken at locations with depth greater than 200 m

setwd('~/Documents/nasa/data/lowres_4km')
write.csv(habs_covar_agg,'habs_covariates_agg.csv',row.names = F)

test <- na.omit(habs_covar_agg)
dim(test)==dim(habs_covar_agg) ### should be true

par(mfrow=c(2,1))
plot(habs_covar_agg$date,habs_covar_agg$CELLCOUNT+1,log='y')
plot(habs$date,habs$CELLCOUNT+1,log='y')


### initial model check
covar <- names(habs_covar_agg)[c(9:22,24:25,27)]

for(i in 1:length(covar)){
  x_y <- formula(paste0('pa100k~',covar[i]))
  mod1 <- glm(x_y,data=habs_covar_agg,family=binomial(link='logit'))
  print(x_y)
  print(summary(mod1))
  cat('\n\n************ANOVA************\n\n')
  print(anova(mod1,test='Chisq'))
  cat('\n\n************END************\n\n')
}

mod1 <- glm(pa100k~as.factor(month)+chlor_a+chl_anom+rbd+nflh+nflh_anom+ssnlw488+cm_bbp+bbp_morel+bbp_carder+abi+rrs_667,
            data=habs_covar_agg,family=binomial(link='logit'))
summary(mod1)
anova(mod1,test='Chisq')

preds <- predict(mod1,newdata=habs_covar_agg,se.fit=T,type='response')
preds$month <- aggregate(preds$fit,by=list(habs_covar_agg$month),mean,na.rm=T)
preds$month.se <- aggregate(preds$se.fit,by=list(habs_covar_agg$month),mean,na.rm=T)
preds$year <- aggregate(preds$fit,by=list(habs_covar_agg$year),mean,na.rm=T)
preds$year.se <- aggregate(preds$se.fit,by=list(habs_covar_agg$year),mean,na.rm=T)

plot(preds$month$Group.1,preds$month$x,pch=18,ylim=c(0,.2))
arrows(preds$month$Group.1,preds$month$x+preds$month.se$x,
         preds$month$Group.1,preds$month$x-preds$month.se$x,length=.105,code=3,angle=90)

plot(preds$year$Group.1,preds$year$x,pch=18,ylim=c(0,.2),typ='n')
# arrows(preds$year$Group.1,preds$year$x+preds$year.se$x,
       # preds$year$Group.1,preds$year$x-preds$year.se$x,length=.105,code=3,angle=90)
polygon(c(preds$year$Group.1,rev(preds$year$Group.1)),
        c(preds$year$x+preds$year.se$x,rev(preds$year$x-preds$year.se$x)),col='gray80')
points(preds$year$Group.1,preds$year$x,pch=18,typ='l')

library(pROC)

temproc <- roc(habs_covar_agg$pa100k , preds$fit, plot=TRUE, grid=TRUE)
# CALCULATE AREA UNDER THE CURVE
temproc$auc  
# Area under the curve: 0.6629
# CONSTRUCT MATRIX OF ROC INFORMATION FOR EACH CUTOFF ("thresholds")	 
roctable <- cbind(temproc$sensitivities, temproc$specificities, temproc$thresholds, 
                  temproc$sensitivities+temproc$specificities)
# FIND CUTOFF WHERE SUM OF THE SENSITIVITY AND SPECIFITY IS MAX
# Sensitivity = proportion of actual positives which are correctly identified as such
# Specificity = proportion of negatives which are correctly identified as such
max(roctable[,1]+roctable[,2])
# [1] 1.241529
# PRINT RECORD FOR THE MAX VALUE TO FIND CUTOFF (= 0.0862855  HERE)
Threshold=roctable[roctable[,4] == max(roctable[,4]),][3]
Threshold 
# [1] 0.07898927
TT=table(mod1$fitted>Threshold, habs_covar_agg$pa100k)
#          0     1
# FALSE 13108   621
# TRUE  11104  1450
FPR =  TT[2,1]/sum(TT[ ,1 ]) 
FNR =   TT[1,2]/sum(TT[ ,2 ])    
FPR # 0.4586156
FNR # 0.2998551

yr <- 2005
subset <- habs_covar_agg[which(habs_covar_agg$year==yr ),]
phat1 <- preds$fit[which(habs_covar_agg$year==yr )]

par(mar=c(5,5,1,6))
plot(subset$LONGITUDE,subset$LATITUDE,asp=1)
quilt.plot(subset$LONGITUDE,subset$LATITUDE,phat1,col=plasma(60),asp=1,add=T)


library(mgcv)

AllModel  <- gam(as.factor(pa100k) ~ as.factor(month) + te(chl_anom,week) + te(chlor_a,week) + te(bbp_carder,bbp_morel) + 
                   te(sst, depth_m) + te(rrs_667,week) + te(LONGITUDE,LATITUDE) + te(rbd,week) + te(depth_m), 
                 data=habs_covar_agg, family = binomial, select=TRUE, method="REML")
save(AllModel, file = "AllModel_initial.RData")
setwd('~/Documents/nasa/data/lowres_4km')
load('AllModel_initial.RData')
p <- plot(AllModel, pages=1, se=TRUE, cex.axis=2, cex.lab=1.5)
summary(AllModel)


