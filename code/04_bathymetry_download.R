library(fields)
library(lubridate)
library(ncdf4)
library(ncdf4.helpers)
library(NISTunits)

vec_brk <- function(input) {
  core <- (diff(input)/2) + input[1:(length(input)-1)]
  beg <- input[1] - (diff(input[1:2])/2)
  last <- input[length(input)] + (diff(input[(length(input)-1):(length(input))])/2)
  output <- c(beg,core,last)
  output
}

### 3) download bathymetry

# a) what bounding box
# (-87.5 30.7, -81, 24.2)
lonbox_w <- -87.5 ### mouth of Mississippi River
latbox_n <- 30.7 ### northern coast
lonbox_e <- -81 ### Florida Bay
latbox_s <- 24.2 ### southern edge of Key West
latitude = c(latbox_s, latbox_n)
longitude = c(lonbox_w, lonbox_e)


# get MODIS lonlat grid
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

### bathymetry
# url <- 'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/60s_bed_elev_netcdf/ETOPO_2022_v1_60s_N90W180_bed.nc'
## higher resolution
url <- 'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/30s/30s_bed_elev_netcdf/ETOPO_2022_v1_30s_N90W180_bed.nc'
data <- nc_open(url)
lon <- ncvar_get(data,'lon')
lon_ind <- which(lon>=lonbox_w & lon<=lonbox_e)
lon_tmp <- lon[lon_ind]
lat <- ncvar_get(data,'lat')
lat_ind <- which(lat>=latbox_s & lat<=latbox_n)
lat_tmp <- lat[lat_ind]
lon_start <- lon_ind[1]
lon_count <- length(lon_ind)
lat_start <- lat_ind[1]
lat_count <- length(lat_ind)
lonlat <- expand.grid(lon=lon_tmp,lat=lat_tmp)
lon_c <- cut(lonlat$lon,vec_brk(lon2))
lat_c <- cut(lonlat$lat,vec_brk(lat2))
lonlat <- expand.grid(lon=levels(lon_c),lat=levels(lat_c))

bathy <- ncvar_get(data,'z',
                   start=c(lon_start,lat_start),
                   count=c(lon_count,lat_count))
bathy[which(bathy>=0)] <- NA
bathy <- (-bathy)
imagePlot(log10(bathy),asp=1,col=rev(viridis(60)),nlevel=59)

bathy_agg <- aggregate(as.vector(bathy),by=list(lon=lon_c,lat=lat_c),mean,na.rm=T)
bathy_agg <- merge(lonlat,bathy_agg,by=c('lon','lat'),all=T)
bathy_agg_m <- t(matrix(bathy_agg$x,length(levels(lat_c)),length(levels(lon_c))))
imagePlot(log10(bathy_agg_m),asp=1,col=rev(viridis(60)),nlevel=59)

hist((bathy_agg_m))
range(bathy_agg_m,na.rm=T)
min(bathy_agg_m,na.rm=T)
max(bathy_agg_m,na.rm=T)
hist((bathy))
range(bathy,na.rm=T)
min(bathy,na.rm=T)
max(bathy,na.rm=T)

setwd('~/Documents/nasa/data/lowres_4km')
saveRDS(bathy_agg_m,'gom_wfs_bathy_lowres.rds') # netcdf is smaller and contains metadata
saveRDS(bathy,'gom_wfs_bathy_hires.rds') # netcdf is smaller and contains metadata
