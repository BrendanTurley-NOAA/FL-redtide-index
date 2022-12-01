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

### average duplicates at same time and location; this disregards ancillary data like temperature, salinity, and wind
habs_sub3 <- habs[,c(3:4,26)]
names(habs)[c(3:4,26)]
ind3 <- duplicated(habs_sub3)
ind3.1 <- duplicated(habs_sub3,fromLast = T)
length(which(ind3))

dups1 <- which(ind3)
dups2 <- which(ind3.1)

plot(habs$date[dups1],habs$CELLCOUNT[dups1],col=2)
points(habs$date[dups2],habs$CELLCOUNT[dups2],col=3)

new <- rep(NA,length(dups1))
for(i in 1:length(dups1)){
  new[i] <- mean(habs$CELLCOUNT[c(dups1[i],dups2[i])],na.rm=T)
}

habs$CELLCOUNT[dups1] <- new
plot(habs$date[dups1],habs$CELLCOUNT[dups1],col=2)
habs <- habs[-dups2,]

### plot data to see if it looks right
plot(habs$LONGITUDE,habs$LATITUDE,asp=1)
plot(habs$date,habs$CELLCOUNT)
plot(habs$date,habs$CELLCOUNT+1,log='y')

habs <- habs[order(habs$date),]
write.csv(habs,'habsos_subset_FL03-21.csv',quote=T,row.names=F)
