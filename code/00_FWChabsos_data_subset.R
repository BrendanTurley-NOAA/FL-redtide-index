library(lubridate)

setwd('~/Desktop/professional/projects/Postdoc_FL/data/habs')

# original data requested from https://habsos.noaa.gov/about
habs1 <- read.csv('habsos_20220225.csv')

habs1$date <- dmy_hm(paste(substr(habs1$SAMPLE_DATE,1,15),substr(habs1$SAMPLE_DATE,30,31)))

# habs <- habs1[which(year(habs1$date)<2023 & year(habs1$date)>1997 & habs1$STATE_ID=='FL' ),]

# plot(habs$LONGITUDE,habs$LATITUDE,asp=1)
# plot(habs$date,habs$CELLCOUNT)
# plot(habs$date,habs$CELLCOUNT+1,log='y')

# write.csv(habs,'habsos_subset_FL98-22.csv')


### updated for matching up with MODIS data
habs <- habs1[which(year(habs1$date)<2022 & year(habs1$date)>2002 & habs1$STATE_ID=='FL' ),]
habs <- habs[which(habs$SAMPLE_DEPTH<=1.5),]

### 4km Aqua MODIS data bounding box
# (-87.5 30.7, -81, 24.2)
lonbox_w <- -87.5 ### mouth of Mississippi River
latbox_n <- 30.7 ### northern coast
lonbox_e <- -81 ### Florida Bay
latbox_s <- 24.2 ### southern edge of Key West
habs <- habs[which(habs$LONGITUDE<=lonbox_e & habs$LONGITUDE>=lonbox_w &
        habs$LATITUDE>=latbox_s & habs$LATITUDE<=latbox_n),]

### exclusions - east coast of FL
ind_ex <- which(habs$LONGITUDE>=(-82) & habs$LATITUDE>=28)
habs <- habs[-ind_ex,]

plot(habs$LONGITUDE,habs$LATITUDE,asp=1)
plot(habs$date,habs$CELLCOUNT)
plot(habs$date,habs$CELLCOUNT+1,log='y')

write.csv(habs,'habsos_subset_FL03-21.csv',quote=F,row.names=F)
