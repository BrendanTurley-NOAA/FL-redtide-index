library(lubridate)
library(scales)

setwd('~/Desktop/professional/projects/Postdoc_FL/data/habs')

# original data requested from https://habsos.noaa.gov/about
habs1 <- read.csv('habsos_20220225.csv')
habs1$date_utc <- dmy_hm(paste(substr(habs1$SAMPLE_DATE,1,15),substr(habs1$SAMPLE_DATE,30,31)))
habs1 <- habs1[which(year(habs1$date_utc)<2021 & year(habs1$date_utc)>1997 & habs1$STATE_ID=='FL' ),]
habs1$dataset <- 'HABSOS'

# original data requested from https://myfwc.com/research/redtide/monitoring/database/
habs2 <- read.csv('FL_Kbrevis_1953-2021.03.08.csv')
names(habs2) <- toupper(names(habs2))
habs2$date <- mdy(habs2$SAMPLE.DATE)
habs2 <- habs2[which(year(habs2$date)<2021 & year(habs2$date)>1997),]
# table(habs2$TIME.ZONE,habs2$COUNTY)
### standardize timezones
switch_data <- function(x){
  renamed <- rep(NA,length(x))
  for(i in 1:length(x)){
    if(nchar(x[i])>1){
      renamed[i]<-  switch(as.character(x[i]),
                           'CDT' = 'America/Chicago',
                           'CST' = 'America/Chicago',
                           'EDT' = 'America/New_York',
                           'EST' = 'America/New_York',
                           'GMT' = 'UTC')
    }
  }
  return(renamed)
}
timez <- switch_data(habs2$TIME.ZONE)
# length(which(is.na(timez)))
### which counties are central timezone?
flco_cen <- c('Bay','Gulf','Escambia','Okaloosa','Santa Rosa','Walton')
### assume no reported time zone is local standard
timez[which(is.na(timez) & is.element(habs2$COUNTY,flco_cen))] <- 'America/Chicago'
timez[which(is.na(timez) & !is.element(habs2$COUNTY,flco_cen))] <- 'America/New_York'
# length(which(is.na(timez)))
# table(timez)
### missing times zones = UTC; assumed
time_missing <- which(nchar(habs2$TIME.ZONE)==0)
habs2$TIME.ZONE[time_missing] <- 'UTC'
timez[time_missing] <- 'UTC'
### missing times = midnight for convenience
time_missing <- which(nchar(habs2$SAMPLE.TIME)==0)
habs2$SAMPLE.TIME[time_missing] <- '00:00'
times <- mdy_hm(paste(habs2$SAMPLE.DATE,habs2$SAMPLE.TIME)) ### assumes everything is UTC and won't accept timezone as a vector
# tz(times)
### apply local timezone, then convert to UTC
habs2$date_utc <- force_tzs(times,timez)
habs2$dataset <- 'FWC'
### standardize depth column name for merge
names(habs2)[grep('depth',names(habs2),ignore.case = T)] <- names(habs1)[grep('depth',names(habs1),ignore.case = T)]


### merge the data on date, lat/lon, and depth
### merge issue so far was due to decimal place reporting
habs1$LATITUDE <- round(habs1$LATITUDE,5)
habs1$LONGITUDE <- round(habs1$LONGITUDE,5)
habs2$LATITUDE <- round(habs2$LATITUDE,5)
habs2$LONGITUDE <- round(habs2$LONGITUDE,5)
habs1$SAMPLE_DEPTH <- round(habs1$SAMPLE_DEPTH,3)
habs2$SAMPLE_DEPTH <- round(habs2$SAMPLE_DEPTH,3)
habs_merge <- merge(habs1,habs2,by=c('date_utc','LATITUDE','LONGITUDE','SAMPLE_DEPTH'),all=T)
# habs_overlap <- merge(habs1,habs2,by=c('date_utc','LATITUDE','LONGITUDE','SAMPLE_DEPTH'),all=F)
overlap <- which(!is.na(habs_merge$dataset.x) & !is.na(habs_merge$dataset.y))
length(overlap)/nrow(habs_merge)
habs_overlap <- habs_merge[overlap,]
### should be one-to-one
plot(habs_overlap$CELLCOUNT+1,
     habs_overlap$KARENIA.BREVIS.ABUNDANCE..CELLS.L.+1,
     log='xy',pch=16,col=alpha(1,.2))
habs_overlap$cell_diff <- (habs_overlap$CELLCOUNT-habs_overlap$KARENIA.BREVIS.ABUNDANCE..CELLS.L.)
length(which(habs_overlap$cell_diff>0 | habs_overlap$cell_diff<0))/nrow(habs_overlap)
hist(habs_overlap$cell_diff)
# plot(habs_overlap$LONGITUDE,
#      habs_overlap$LATITUDE,
#      asp=1,cex=log(abs(habs_overlap$cell_diff))/5,
#      pch=16,col=alpha(1,.1))
par(mfrow=c(1,2))
plot(habs_overlap$LONGITUDE[which(habs_overlap$cell_diff<0)],
       habs_overlap$LATITUDE[which(habs_overlap$cell_diff<0)],
       cex=log(abs(habs_overlap$cell_diff[which(habs_overlap$cell_diff<0)]))/5,
     asp=1,pch=21,col=alpha(1,.1),bg=alpha('blue',.1))
plot(habs_overlap$LONGITUDE[which(habs_overlap$cell_diff>0)],
       habs_overlap$LATITUDE[which(habs_overlap$cell_diff>0)],
       cex=log(abs(habs_overlap$cell_diff[which(habs_overlap$cell_diff>0)]))/5,
       asp=1,pch=21,col=alpha(1,.1),bg=alpha('red',.1))
### these plots are the same; what is going on?
length(which(habs_overlap$cell_diff>0 | habs_overlap$cell_diff<0))
rbind(habs_overlap[which.max(habs_overlap$cell_diff),],
      habs_overlap[which.min(habs_overlap$cell_diff),])
habs_overlap_diff <- habs_overlap[which(habs_overlap$cell_diff>0 | habs_overlap$cell_diff<0),]
habs_overlap_diff<- habs_overlap_diff[order(habs_overlap_diff$date_utc,habs_overlap_diff$LATITUDE,habs_overlap_diff$LONGITUDE),]
plot(habs_overlap_diff$cell_diff[seq(1,nrow(habs_overlap_diff),2)],
     habs_overlap_diff$cell_diff[seq(2,nrow(habs_overlap_diff),2)])

dups <- which(duplicated(habs_overlap_diff))
which(duplicated(habs_overlap_diff,fromLast = T))
habs_overlap_diff <- habs_overlap_diff[-dups,]
# habs_overlap_diff[127:128,]
dups1 <- which(duplicated(habs_overlap_diff[,1:3]))
dups2 <- which(duplicated(habs_overlap_diff[,1:3],fromLast = T))
reps <- data.frame(dups1, habs_overlap_diff$CELLCOUNT[dups1],
                   dups2, habs_overlap_diff$KARENIA.BREVIS.ABUNDANCE..CELLS.L.[dups2])
names(reps) <- c('inx2','habsos_cells','inx2','fwc_cells')
reps$diff <- reps$habsos_cells-reps$fwc_cells
inx <- which(reps$diff!=0)
habs_dups <- habs_overlap_diff[sort(c(reps[inx,1],reps[inx,3])),]

habsos <- which(habs_merge$dataset.x=='HABSOS')
habsos_r <- which(habs_merge$dataset.x=='HABSOS' & is.na(habs_merge$dataset.y))
fwc <- which(habs_merge$dataset.y=='FWC')
fwc_r <- which(habs_merge$dataset.y=='FWC' & is.na(habs_merge$dataset.x))

length(overlap) + length(habsos_r) + length(fwc_r)

#### investigating merge
fwc_r_comp <- rep(NA,length(habsos_r))
fwc_r_comp[1:length(fwc_r)] <- fwc_r
comps <- data.frame(habsos=habsos_r,fwc=fwc_r_comp)

plot(habsos_r,fwc_r_comp)
plot(habsos_r-fwc_r_comp)

rbind(habs_merge[habsos_r[1],],habs_merge[fwc_r[1],])

habs_merge$LATITUDE[fwc_r[1]]==habs_merge$LATITUDE[habsos_r[1]]
habs_merge$LATITUDE[fwc_r[1]]-habs_merge$LATITUDE[habsos_r[1]] #WTF
habs_merge$LONGITUDE[fwc_r[1]]-habs_merge$LONGITUDE[habsos_r[1]] #WTF
habs_merge$SAMPLE_DEPTH[fwc_r[1]]-habs_merge$SAMPLE_DEPTH[habsos_r[1]] 

# library(dplyr)
# habs_merge2 <- full_join(habs1,habs2,by=c('date_utc','LATITUDE','LONGITUDE','SAMPLE_DEPTH'))
### same result


### lots of overlap not merge because time of sampling not matching
par(mfrow=c(1,2))
plot(habs_merge$LONGITUDE[overlap],habs_merge$LATITUDE[overlap],asp=1)
points(habs_merge$LONGITUDE[habsos_r],habs_merge$LATITUDE[habsos_r],col='red',pch=16,cex=.2)
plot(habs_merge$LONGITUDE[overlap],habs_merge$LATITUDE[overlap],asp=1)
points(habs_merge$LONGITUDE[fwc_r],habs_merge$LATITUDE[fwc_r],col='green',pch=16,cex=.2)




### merge the data on date, lat/lon, and depth
### strip time off
habs1$date_utc[1]
habs1$date_ymd <- ymd(substr(habs1$date_utc,1,10))
habs2$date_ymd <- ymd(substr(habs2$date_utc,1,10))
habs_merge <- merge(habs1,habs2,by=c('date_ymd','LATITUDE','LONGITUDE','SAMPLE_DEPTH'),all=T)
overlap <- which(!is.na(habs_merge$dataset.x) & !is.na(habs_merge$dataset.y))
length(overlap)/nrow(habs_merge)
habs_overlap <- habs_merge[overlap,]
### should be one-to-one
plot(habs_overlap$CELLCOUNT+1,habs_overlap$KARENIA.BREVIS.ABUNDANCE..CELLS.L.+1,log='xy')
habs_overlap$cell_diff <- (habs_overlap$CELLCOUNT-habs_overlap$KARENIA.BREVIS.ABUNDANCE..CELLS.L.)
length(which(habs_overlap$cell_diff>0))/nrow(habs_overlap)
hist(habs_overlap$cell_diff)
plot(habs_overlap$LONGITUDE,habs_overlap$LATITUDE,asp=1,cex=log(abs(habs_overlap$cell_diff))/5)

habsos <- which(habs_merge$dataset.x=='HABSOS')
habsos_r <- which(habs_merge$dataset.x=='HABSOS' & is.na(habs_merge$dataset.y))
fwc <- which(habs_merge$dataset.y=='FWC')
fwc_r <- which(habs_merge$dataset.y=='FWC' & is.na(habs_merge$dataset.x))

length(overlap) + length(habsos_r) + length(fwc_r)

### lots of overlap not merge because time of sampling not matching
par(mfrow=c(1,2))
plot(habs_merge$LONGITUDE[overlap],habs_merge$LATITUDE[overlap],asp=1)
points(habs_merge$LONGITUDE[habsos_r],habs_merge$LATITUDE[habsos_r],col='red',pch='.')
plot(habs_merge$LONGITUDE[overlap],habs_merge$LATITUDE[overlap],asp=1)
points(habs_merge$LONGITUDE[fwc_r],habs_merge$LATITUDE[fwc_r],col='blue',pch='.')
