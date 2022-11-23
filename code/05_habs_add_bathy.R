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

setwd('~/Desktop/professional/projects/Postdoc_FL/data/habs')
habs <- read.csv('habsos_subset_FL03-21.csv')
habs$date <- ymd_hms(habs$date)
habs <- habs[order(habs$date),]

setwd('~/Documents/nasa/data/lowres_4km')
load('bathy.RData')

### cut by bathymetry grid
habs$lon_b <- cut(habs$LONGITUDE,vec_brk(lon_bathy))
habs$lat_b <- cut(habs$LATITUDE,vec_brk(lat_bathy))

bathy <- data.frame(longitude=lonlat_bathy$lon,
                    latitude=lonlat_bathy$lat,
                    depth_m=as.vector(bathy_orig))
bathy$lon_b <- cut(bathy$longitude,vec_brk(lon_bathy))
bathy$lat_b <- cut(bathy$latitude,vec_brk(lat_bathy))

hab_bathy <- merge(habs,bathy[,-c(1,2)],by=c('lon_b','lat_b'),all.x=T) # don't include superfluous lon/lats
write.csv(hab_bathy,'habs_bathy.csv',quote=T,row.names=F)


### examine the results
ind <- which(hab_bathy$depth_m>=1)
which(hab_bathy$depth_m=='NaN')
hist(hab_bathy$SAMPLE_DEPTH[ind])
hist(hab_bathy$depth_m[ind])

plot(hab_bathy$LONGITUDE,hab_bathy$LATITUDE,cex=log10(hab_bathy$depth_m))
points(hab_bathy$LONGITUDE[ind],hab_bathy$LATITUDE[ind],col=2)

plot(hab_bathy$LONGITUDE[ind],hab_bathy$LATITUDE[ind])

plot(hab_bathy$SAMPLE_DEPTH[ind],hab_bathy$depth_m[ind])
plot(hab_bathy$date[ind],hab_bathy$CELLCOUNT[ind])
plot(hab_bathy$date[ind],hab_bathy$CELLCOUNT[ind]+1,log='y')
plot(year(hab_bathy$date[ind]),hab_bathy$CELLCOUNT[ind]+1,log='y')
plot(month(hab_bathy$date[ind]),hab_bathy$CELLCOUNT[ind]+1,log='y')


bathy_0 <- bathy_orig
bathy_0[which(bathy_0<=0)] <- NA
bathy_0[which(bathy_0>5)] <- 5

png('Florida.png',width=10,height=10,units='in',res=500)
imagePlot(lon_bathy,lat_bathy,bathy_0,asp=1)
dev.off()

imagePlot(bathy_agg_modis,asp=1)
bathy_0 <- bathy_agg_modis
bathy_0[which(bathy_0<=0)] <- NA
imagePlot(bathy_0,asp=1)
