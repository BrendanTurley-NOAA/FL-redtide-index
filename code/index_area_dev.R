gc()

library(fields)
library(MASS)
library(lubridate)
library(scales)

setwd('~/Documents/nasa/data/lowres_4km')
# write.csv(habs_covar_agg,'habs_covariates_agg.csv',row.names = F)
habs_covar_agg <- read.csv('habs_covariates_agg.csv')
habs_covar_agg$date <- ymd(habs_covar_agg$date)

habs <- kde2d(habs_covar_agg$LONGITUDE,habs_covar_agg$LATITUDE,n=158)

threshold <- quantile(habs$z,.5)
habs$z[which(habs$z<threshold)] <- NA

imagePlot(habs,asp=1,nlevel=64,col=inferno(63))
points(habs_covar_agg$LONGITUDE,habs_covar_agg$LATITUDE,col=4,pch=20,cex=.5)

hull <- chull(habs_covar_agg$LONGITUDE,habs_covar_agg$LATITUDE)
centroid <- cbind(lon=weighted.mean(habs_covar_agg$LONGITUDE,habs_covar_agg$CELLCOUNT),
                  lat=weighted.mean(habs_covar_agg$LATITUDE,habs_covar_agg$CELLCOUNT))
centroid2 <- cbind(lon=mean(habs_covar_agg$LONGITUDE),
                  lat=mean(habs_covar_agg$LATITUDE))

plot(habs_covar_agg$LONGITUDE,habs_covar_agg$LATITUDE,col=4,pch=20,cex=.5,asp=1)
points(habs_covar_agg$LONGITUDE[c(hull,hull[1])],
       habs_covar_agg$LATITUDE[c(hull,hull[1])],typ='o')
points(centroid,col=2)
points(centroid2,col=3)

lon_dist <- sapply(habs_covar_agg$LONGITUDE,function(x) x-centroid[1])
lat_dist <- sapply(habs_covar_agg$LATITUDE,function(x) x-centroid[2])

dists <- sqrt(lon_dist^2 + lat_dist^2)
threshold <- quantile(dists,.75)

imcp75 <- which(dists<threshold)
mcp75 <- habs_covar_agg[imcp75,]
hull2 <- chull(mcp75$LONGITUDE,mcp75$LATITUDE)

plot(habs_covar_agg$LONGITUDE,habs_covar_agg$LATITUDE,col=4,pch=20,cex=.5,asp=1)
points(habs_covar_agg$LONGITUDE[c(hull,hull[1])],
       habs_covar_agg$LATITUDE[c(hull,hull[1])],typ='o',col=4)
points(mcp75$LONGITUDE,mcp75$LATITUDE,col=2,pch=20,cex=.5)
points(mcp75$LONGITUDE[c(hull2,hull2[1])],
       mcp75$LATITUDE[c(hull2,hull2[1])],typ='o',col=2)
points(centroid,col=3)

wss <- rep(NA,10)
for(i in 1:10){
  km <- kmeans(habs_covar_agg[,c(5,4)],i,nstart=4)
  wss[i] <- km$tot.withinss
}
plot(wss,typ='b')

k <- 6
km <- kmeans(habs_covar_agg[,c(5,4)],k,nstart=4)

plot(habs_covar_agg$LONGITUDE,habs_covar_agg$LATITUDE,col='gray50',pch=20,cex=.5,asp=1)
# points(km$centers)
# points(habs_covar_agg$LONGITUDE,habs_covar_agg$LATITUDE,col=km$cluster,pch=20,cex=.5)

hulls <- matrix(NA,50*k,3)
m <- 1
n <- 0
for(i in 1:k){
  ind <- which(km$cluster==i)
  tmp <- habs_covar_agg[ind,]
  
  centroid <- cbind(lon=weighted.mean(tmp$LONGITUDE,tmp$CELLCOUNT),
                    lat=weighted.mean(tmp$LATITUDE,tmp$CELLCOUNT))
  
  lon_dist <- sapply(tmp$LONGITUDE,function(x) x-centroid[1])
  lat_dist <- sapply(tmp$LATITUDE,function(x) x-centroid[2])
  
  dists <- sqrt(lon_dist^2 + lat_dist^2)
  threshold <- quantile(dists,.95)
  
  imcp75 <- which(dists<threshold)
  mcp75 <- tmp[imcp75,]
  hull <- chull(mcp75$LONGITUDE,mcp75$LATITUDE)
  n <- n + length(hull)
  hulls[m:n,] <- cbind(i,mcp75$LONGITUDE[hull],mcp75$LATITUDE[hull])
  m <- n + 1
  
  points(mcp75$LONGITUDE[c(hull,hull[1])],mcp75$LATITUDE[c(hull,hull[1])],col=i,typ='o')
}

