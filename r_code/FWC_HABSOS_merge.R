library(lubridate)
library(scales)

setwd('~/Desktop/professional/projects/Postdoc_FL/data/habs')

# original data requested from https://habsos.noaa.gov/about
habs1 <- read.csv('habsos_20220225.csv')
habs1$date_utc <- dmy_hm(paste(substr(habs1$SAMPLE_DATE,1,15),substr(habs1$SAMPLE_DATE,30,31)))
habs1 <- habs1[which(year(habs1$date_utc)<2021 & year(habs1$date_utc)>1997 & habs1$STATE_ID=='FL' ),]
habs1$dataset <- 'HABSOS'
### duplicates?
# dups <- which(duplicated(habs1))
# dupr <- which(duplicated(habs1,fromLast = T))

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
# habs2 <- habs2[-time_missing,]
# timez <- timez[-time_missing]
habs2$SAMPLE.TIME[time_missing] <- '00:00'
times <- mdy_hm(paste(habs2$SAMPLE.DATE,habs2$SAMPLE.TIME)) ### assumes everything is UTC and won't accept timezone as a vector
# tz(times)
### apply local timezone, then convert to UTC
habs2$date_utc <- force_tzs(times,timez)
habs2$dataset <- 'FWC'
### duplicates?
dups <- which(duplicated(habs2))
# dupr <- which(duplicated(habs2,fromLast = T))
habs2 <- habs2[-dups,]


### merge the data on date, lat/lon, and depth
### standardize depth column name for merge
names(habs2)[grep('depth',names(habs2),ignore.case = T)] <- names(habs1)[grep('depth',names(habs1),ignore.case = T)]
### merge issue so far was due to decimal place reporting
habs1$LATITUDE <- round(habs1$LATITUDE,5)
habs1$LONGITUDE <- round(habs1$LONGITUDE,5)
habs2$LATITUDE <- round(habs2$LATITUDE,5)
habs2$LONGITUDE <- round(habs2$LONGITUDE,5)
habs1$SAMPLE_DEPTH <- round(habs1$SAMPLE_DEPTH,3)
habs2$SAMPLE_DEPTH <- round(habs2$SAMPLE_DEPTH,3)
habs_merge <- merge(habs1,habs2,by=c('date_utc','LATITUDE','LONGITUDE','SAMPLE_DEPTH'),all=T)
### duplicates?
# dups <- which(duplicated(habs_merge))
# dupr <- which(duplicated(habs_merge,fromLast = T))
### what is the overlap?
# habs_overlap <- merge(habs1,habs2,by=c('date_utc','LATITUDE','LONGITUDE','SAMPLE_DEPTH'),all=F)
overlap <- which(!is.na(habs_merge$dataset.x) & !is.na(habs_merge$dataset.y))
length(overlap)/nrow(habs_merge)
habs_overlap <- habs_merge[overlap,]
### duplicates?
# dups <- which(duplicated(habs_merge))
# dupr <- which(duplicated(habs_merge,fromLast = T))
inx <- which(duplicated(habs_overlap$OBJECTID))
iny <- which(duplicated(habs_overlap$OBJECTID,fromLast = T))
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

### these plots are bascially the same; what is going on?
### many samples at same time and location, but merge fxn can't discriminate
length(which(habs_overlap$cell_diff>0 | habs_overlap$cell_diff<0))
rbind(habs_overlap[which.max(habs_overlap$cell_diff),],
      habs_overlap[which.min(habs_overlap$cell_diff),])

habs_overlap_diff <- habs_overlap[which(habs_overlap$cell_diff>0 | habs_overlap$cell_diff<0),]
habs_overlap_diff<- habs_overlap_diff[order(habs_overlap_diff$date_utc,habs_overlap_diff$LATITUDE,habs_overlap_diff$LONGITUDE),]
habs_overlap_diff$dif2 <- abs(habs_overlap_diff$cell_diff)
dups <- which(duplicated(habs_overlap_diff[,c(1:4,46)]))
dupr <- which(duplicated(habs_overlap_diff[,c(1:4,46)],fromLast = T))
plot(dups-dupr)
plot(diff(sort(c(dups,dupr))))

### how many samples per unique timedate, location, sample depth?
test_dup <- paste(habs_overlap_diff$date_utc,habs_overlap_diff$LATITUDE,habs_overlap_diff$LONGITUDE,habs_overlap_diff$SAMPLE_DEPTH,sep=' ')
sort(table(test_dup))
habs_overlap_diff[which(habs_overlap_diff$date_utc==ymd_hms('2003-01-17 12:00:00')),]


dups1 <- which(duplicated(habs_overlap_diff[,1:3]))
dups2 <- which(duplicated(habs_overlap_diff[,1:3],fromLast = T))
plot(dups1-dups2)
plot(diff(sort(c(dups1,dups2))))
reps <- data.frame(dups1, habs_overlap_diff$CELLCOUNT[dups1],
                   dups2, habs_overlap_diff$KARENIA.BREVIS.ABUNDANCE..CELLS.L.[dups2])
names(reps) <- c('inx1','habsos_cells','inx2','fwc_cells')
inx_dups1 <- which(is.element(reps$inx1,reps$inx2))
inx_dups2 <- which(is.element(reps$inx2,reps$inx1))
habs_overlap_diff[103:107,]

reps$diff <- reps$habsos_cells-reps$fwc_cells
inx <- which(reps$diff!=0)
habs_dups <- habs_overlap_diff[sort(c(reps[inx,1],reps[inx,3])),]


### separate non-overlapping
habsos <- which(habs_merge$dataset.x=='HABSOS')
habsos_r <- which(habs_merge$dataset.x=='HABSOS' & is.na(habs_merge$dataset.y))
fwc <- which(habs_merge$dataset.y=='FWC')
fwc_r <- which(habs_merge$dataset.y=='FWC' & is.na(habs_merge$dataset.x))
### should be TRUE
(length(overlap) + length(habsos_r) + length(fwc_r)) == nrow(habs_merge)

### plot out to see overlap
par(mfrow=c(1,2))
plot(habs_merge$LONGITUDE[overlap],habs_merge$LATITUDE[overlap],asp=1)
points(habs_merge$LONGITUDE[habsos_r],habs_merge$LATITUDE[habsos_r],col='red',pch=16,cex=.2)
plot(habs_merge$LONGITUDE[overlap],habs_merge$LATITUDE[overlap],asp=1)
points(habs_merge$LONGITUDE[fwc_r],habs_merge$LATITUDE[fwc_r],col='green',pch=16,cex=.2)


habs_combined <- habs_merge[habsos,]
habs_combined2 <- merge(habs1,habs2,by=c('date_utc','LATITUDE','LONGITUDE','SAMPLE_DEPTH'),all.x=T)
dups <- which(duplicated(habs_combined$OBJECTID))
dupr <- which(duplicated(habs_combined$OBJECTID,fromLast = T))
plot(dups-dupr)
habs_combined[2134:2135,]

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



### ECOHAB cruises
ec_in <- grep('ecohab',habs2$COLLECTION.AGENCY,ignore.case = T)
ecohab <- habs2[ec_in,]

plot(ecohab$LONGITUDE,ecohab$LATITUDE,asp=1,col=year(ecohab$date_utc))
plot(ecohab$date_utc,ecohab$KARENIA.BREVIS.ABUNDANCE..CELLS.L.)

table(year(ecohab$date_utc),month(ecohab$date_utc))

cuts <- cut(tmp$KARENIA.BREVIS.ABUNDANCE..CELLS.L.,breaks=c(0,1e3,1e4,1e5,1e6,1e10))

yrs <- sort(unique(year(ecohab$date_utc)))
for(i in yrs){
  tmp <- ecohab[which(year(ecohab$date_utc)==i),]
  tmp <- tmp[order(tmp$KARENIA.BREVIS.ABUNDANCE..CELLS.L.,decreasing = F),]
  cuts <- cut(tmp$KARENIA.BREVIS.ABUNDANCE..CELLS.L.,breaks=c(0,1e3,1e4,1e5,1e6,1e10))
  plot(ecohab$LONGITUDE,ecohab$LATITUDE,asp=1,col='white')
  points(tmp$LONGITUDE,tmp$LATITUDE,asp=1,pch=21,
         # bg=alpha('red',log(tmp$KARENIA.BREVIS.ABUNDANCE..CELLS.L.)/max(log(tmp$KARENIA.BREVIS.ABUNDANCE..CELLS.L.))))
  bg=alpha(c('gray50','white','gold','orange2','red3')[cuts],.7))
  
  mtext(i)
}
