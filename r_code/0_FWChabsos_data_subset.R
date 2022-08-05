library(lubridate)

setwd('~/Desktop/professional/projects/Postdoc_FL/data/habs')

# original data requested from https://habsos.noaa.gov/about
habs1 <- read.csv('habsos_20220225.csv')
# habs2 <- read.csv('habsos_proofed_20211210.csv') # unused; older dataset

habs1$date <- dmy_hm(paste(substr(habs1$SAMPLE_DATE,1,15),substr(habs1$SAMPLE_DATE,30,31)))

habs <- habs1[which(year(habs1$date)<2023 & year(habs1$date)>1997 & habs1$STATE_ID=='FL' ),]

plot(habs$LONGITUDE,habs$LATITUDE,asp=1)
plot(habs$date,habs$CELLCOUNT+1,log='y')

write.csv(habs,'habsos_subset_FL98-22.csv')

