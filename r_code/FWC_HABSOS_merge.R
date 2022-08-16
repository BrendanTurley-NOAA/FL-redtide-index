library(lubridate)

setwd('~/Desktop/professional/projects/Postdoc_FL/data/habs')

# original data requested from https://habsos.noaa.gov/about
habs1 <- read.csv('habsos_20220225.csv')
habs1$date_utc <- dmy_hm(paste(substr(habs1$SAMPLE_DATE,1,15),substr(habs1$SAMPLE_DATE,30,31)))
habs1 <- habs1[which(year(habs1$date_utc)<2023 & year(habs1$date_utc)>1997 & habs1$STATE_ID=='FL' ),]
habs1$dataset <- 'HABSOS'

# original data requested from https://myfwc.com/research/redtide/monitoring/database/
habs2 <- read.csv('FL_Kbrevis_1953-2021.03.08.csv')
names(habs2) <- toupper(names(habs2))
habs2$date <- mdy(habs2$SAMPLE.DATE)
habs2 <- habs2[which(year(habs2$date)<2023 & year(habs2$date)>1997),]
table(habs2$TIME.ZONE,habs2$COUNTY)
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
### missing times = noon for convenience
time_missing <- which(nchar(habs2$SAMPLE.TIME)==0)
habs2$SAMPLE.TIME[time_missing] <- '12:00'
times <- mdy_hm(paste(habs2$SAMPLE.DATE,habs2$SAMPLE.TIME)) ### assumes everything is UTC and won't accept timezone as a vector
# tz(times)
### apply local timezone, then conver to UTC
habs2$date_utc <- force_tzs(times,timez)
habs2$dataset <- 'FWC'
names(habs2)[grep('depth',names(habs2),ignore.case = T)] <- names(habs1)[grep('depth',names(habs1),ignore.case = T)]


habs_merge <- merge(habs1,habs2,by=c('date_utc','LATITUDE','LONGITUDE','SAMPLE_DEPTH'),all=T)
overlap <- which(!is.na(habs_merge$dataset.x) & !is.na(habs_merge$dataset.y))
length(overlap)/nrow(habs_merge)
habs_overlap <- habs_merge[overlap,]
### should be one-to-one
plot(habs_overlap$CELLCOUNT,habs_overlap$KARENIA.BREVIS.ABUNDANCE..CELLS.L.)
habs_overlap$cell_diff <- (habs_overlap$CELLCOUNT-habs_overlap$KARENIA.BREVIS.ABUNDANCE..CELLS.L.)
length(which(habs_overlap$cell_diff>0))/nrow(habs_overlap)
hist(habs_overlap$cell_diff)
plot(habs_overlap$LONGITUDE,habs_overlap$LATITUDE,asp=1,cex=log(abs(habs_overlap$cell_diff))/5)

habsos <- which(habs_merge$dataset.x=='HABSOS')
fwc <- which(habs_merge$dataset.y=='FWC')

plot(habs_merge$LONGITUDE[overlap],habs_merge$LATITUDE[overlap],asp=1)
points(habs_merge$LONGITUDE[habsos],habs_merge$LATITUDE[habsos],col='red',pch='.')
points(habs_merge$LONGITUDE[fwc],habs_merge$LATITUDE[fwc],col='blue',pch='.')
