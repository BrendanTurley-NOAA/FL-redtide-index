library(sp)
library(spatstat)

?point.in.polygon
?ppp
?tess

A <- tess(xgrid=0:4,ygrid=0:4)
A
plot(A)
B <- A[c(1, 2, 5, 7, 9)]
B
plot(B)
v <- as.im(function(x,y){factor(round(5 * (x^2 + y^2)))}, W=owin())
levels(v) <- letters[seq(length(levels(v)))]
E <- tess(image=v)
plot(E)
G <- tess(image=v, marks=toupper(levels(v)), unitname="km")
G
plot(G)

coords <- data.frame(x=sample(seq(0,4,.1),100,replace=T),y=sample(seq(0,4,.1),100,replace=T))
z <- rlnorm(100,0,1)
coords2 <- ppp(coords$x,coords$y,c(0,4),c(0,4))
coords3 <- cut(coords2,A)

plot(coords)
plot(coords2)
plot(coords3)
plot(A,add=T)

zz <- aggregate(z,by=list(coords3$marks),mean,na.rm=T)

lonlat <- tess(xgrid=lon2,ygrid=sort(lat2))

cord.dec = SpatialPoints(cbind(lon2, lat2), proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS('+proj=utm +zone=17 +units=km'))
plot(cord.UTM,axes=T)



url <- 'http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2011/002/A2011002.L3m_DAY_RRS_Rrs_443_4km.nc'
# system.time(modis1 <- nc_open(url,readunlim=T,suppress_dimvals=F,return_on_error=T))
system.time(modis1 <- nc_open(url,readunlim=F,suppress_dimvals=T,return_on_error=T))
atts <- ncatt_get(modis1,0)
# names(atts)[-c(1,6,8,9,11:14,18:26,29:34,36:44,46:48,51,53,55,57:64)]
lon <- ncvar_get(modis1, 'lon')
lon_start <- which(lon>lonbox_w)[1]-1
lon_stop <- which(lon>lonbox_e)[1]
lon_count <- length(lon_start:lon_stop)
lat <- ncvar_get(modis1, 'lat')
lat_start <- which(lat<latbox_n)[1]-1
lat_stop <- which(lat<latbox_s)[1]
lat_count <- length(lat_start:lat_stop)
lon2 <- ncvar_get(modis1, 'lon',start=lon_start,count=lon_count)
lat2 <- ncvar_get(modis1, 'lat',start=lat_start,count=lat_count)

lonlat <- expand.grid(lon=lon2,lat=lat2)
# lonlat <- expand.grid(lon=lon2[sample(1:length(lon2),10)],lat=lat2[sample(1:length(lat2),10)])
z <- rlnorm(nrow(lonlat),0,1)

system.time({
  lon_c <- cut(lonlat$lon,vec_brk(lon2))
  lat_c <- cut(lonlat$lat,vec_brk(lat2))
  lonlat_c <- expand.grid(lon=levels(lon_c),lat=levels(lat_c))
  # z <- rlnorm(nrow(lonlat_c),0,1)
  zz1 <- aggregate(z,by=list(lon=lon_c,lat=lat_c),mean,na.rm=T)
})

system.time({
  lonlat_g <- tess(xgrid=vec_brk(lon2),ygrid=sort(vec_brk(lat2)))
  lonlat_g$marks <- 1:lonlat_g$n
  lonlat_p <- ppp(lonlat$lon,lonlat$lat,range(vec_brk(lon2)),range(vec_brk(lat2)))
  lonlat_pg <- cut(lonlat_p,lonlat_g)
  lonlat_pg$marks <- 1:lonlat_pg$n
  levels(lonlat_pg$marks) <- 1:length(levels(lonlat_pg$marks))
  # z <- rlnorm(lonlat_pg$n,0,1)
  zz2 <- aggregate(z,by=list(lonlat_pg$marks),mean,na.rm=T)
})

x1 <- matrix(zz1$x,158,158)
x2 <- matrix(zz2$x,158,158)

par(mfrow=c(1,2))
imagePlot(x1[,ncol(x1):1],asp=1,col=inferno(60),nlevel=59)
mtext('zz1')
imagePlot(x2,asp=1,col=inferno(60),nlevel=59)
mtext('zz2')

imagePlot(x1[,ncol(x1):1]-x2)
range((x1[,ncol(x1):1]-x2))





plot(lonlat_p)
plot(lonlat_p,xlim=c(-82.5,-82),ylim=c(25,25.5))
plot(lonlat_g,add=T)

plot(lonlat)
plot(lonlat_p)
