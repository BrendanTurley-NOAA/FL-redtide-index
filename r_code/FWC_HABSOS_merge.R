library(lubridate)

setwd('~/Desktop/professional/projects/Postdoc_FL/data/habs')

# original data requested from https://habsos.noaa.gov/about
habs1 <- read.csv('habsos_20220225.csv')
# original data requested from https://myfwc.com/research/redtide/monitoring/database/
habs2 <- read.csv('FL_Kbrevis_1953-2021.03.08.csv')

### issue: not clear what the timezones are for the HABSOS data and what to do about the missing timezones in the FWRI data

habs1$date <- dmy_hm(paste(substr(habs1$SAMPLE_DATE,1,15),substr(habs1$SAMPLE_DATE,30,31)))


table(habs2$Time.Zone)

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
timez[is.na(timez)] <- 'GMT' # placeholder
table(timez)

times <- mdy_hm(paste(habs2$Sample.Date,habs2$Sample.Time,timez))
tz(times)
table(timez[is.na(times)]) # which failed
