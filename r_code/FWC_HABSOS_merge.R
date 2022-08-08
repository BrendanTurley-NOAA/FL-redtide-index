library(lubridate)

setwd('~/Desktop/professional/projects/Postdoc_FL/data/habs')

# original data requested from https://habsos.noaa.gov/about
habs1 <- read.csv('habsos_20220225.csv')
# original data requested from https://myfwc.com/research/redtide/monitoring/database/
habs2 <- read.csv('FL_Kbrevis_1953-2021.03.08.csv')
habs2$date <- mdy(habs2$Sample.Date)
habs2 <- habs2[which(year(habs2$date)<2023 & year(habs2$date)>1997),]

### issue: not clear what the timezones are for the HABSOS data and what to do about the missing timezones in the FWRI data

habs1$date <- dmy_hm(paste(substr(habs1$SAMPLE_DATE,1,15),substr(habs1$SAMPLE_DATE,30,31)))


table(habs2$Time.Zone,habs2$County)


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

timez <- switch_data(habs2$Time.Zone)
length(is.na(timez))
### which counties are central timezone?
flco_cen <- c('Bay','Gulf','Escambia','Okaloosa','Santa Rosa','Walton')
### assume no reported time zone is local standard
timez[which(is.na(timez) & is.element(habs2$County,flco_cen))] <- 'America/Chicago'
timez[which(is.na(timez) & !is.element(habs2$County,flco_cen))] <- 'America/New_York'

table(timez)

times <- mdy_hm(paste(habs2$Sample.Date,habs2$Sample.Time))
tz(times)
table(timez[is.na(times)]) # which failed
