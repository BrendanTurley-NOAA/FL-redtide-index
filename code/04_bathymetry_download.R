library(fields)
library(lubridate)
library(ncdf4)
library(ncdf4.helpers)
library(NISTunits)

lonlat_ind <- function(lonlat,left,right){
  ind <- which(lonlat>=left & lonlat<=right)
  out <- lonlat[ind]
  
  if(out[1]>left){
    ind <- c((ind[1]-1),ind)
    out <- lonlat[ind]
  }
  
  if(out[length(ind)]<right){
    ind <- c(ind,(ind[length(out)]+1))
    out <- lonlat[ind]
  }
  return(list(data=out,indices=ind))
}

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
url <- 'http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2011/0102/AQUA_MODIS.20110102.L3m.DAY.RRS.Rrs_443.4km.nc'
# url <- 'http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2011/002/A2011002.L3m_DAY_RRS_Rrs_443_4km.nc'
# system.time(modis1 <- nc_open(url,readunlim=T,suppress_dimvals=F,return_on_error=T))
system.time(modis1 <- nc_open(url,readunlim=F,suppress_dimvals=T,return_on_error=T))
# atts <- ncatt_get(modis1,0)
# # names(atts)[-c(1,6,8,9,11:14,18:26,29:34,36:44,46:48,51,53,55,57:64)]
# lon <- ncvar_get(modis1, 'lon')
# lon_start <- which(lon>lonbox_w)[1]-1
# lon_stop <- which(lon>lonbox_e)[1]
# lon_count <- length(lon_start:lon_stop)
# lat <- ncvar_get(modis1, 'lat')
# lat_start <- which(lat<latbox_n)[1]-1
# lat_stop <- which(lat<latbox_s)[1]
# lat_count <- length(lat_start:lat_stop)
# lon_modis <- ncvar_get(modis1, 'lon',start=lon_start,count=lon_count)
# lat_modis <- ncvar_get(modis1, 'lat',start=lat_start,count=lat_count)
lon <- ncvar_get(modis1, 'lon')
lons <- lonlat_ind(lon,lonbox_w,lonbox_e)
lon_modis <- lons$data
lat <- ncvar_get(modis1,'lat')
lats <- lonlat_ind(lat,latbox_s,latbox_n)
lat_modis <- lats$data


### bathymetry
### NOAA ETOPO 2022; https://www.ncei.noaa.gov/products/etopo-global-relief-model
## 30-arc seconds ~1km^2
url <- 'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/30s/30s_bed_elev_netcdf/ETOPO_2022_v1_30s_N90W180_bed.nc'
data <- nc_open(url)
lon <- ncvar_get(data,'lon')
lons <- lonlat_ind(lon,lonbox_w,lonbox_e)
lon_ind <- lons$indices
lon_bathy <- lons$data
# lon_ind <- which(lon>=lonbox_w & lon<=lonbox_e)
# lon_bathy <- lon[lon_ind]
lat <- ncvar_get(data,'lat')
lats <- lonlat_ind(lat,latbox_s,latbox_n)
lat_ind <- lats$indices
lat_bathy <- lats$data
# lat_ind <- which(lat>=latbox_s & lat<=latbox_n)
# lat_bathy <- lat[lat_ind]
lon_start <- lon_ind[1]
lon_count <- length(lon_ind)
lat_start <- lat_ind[1]
lat_count <- length(lat_ind)
lonlat_bathy <- expand.grid(lon=lon_bathy,lat=lat_bathy)
lon_c <- cut(lonlat_bathy$lon,vec_brk(lon_modis))
lat_c <- cut(lonlat_bathy$lat,vec_brk(lat_modis))
lonlat_modis <- expand.grid(lon=levels(lon_c),lat=levels(lat_c))

bathy_orig <- ncvar_get(data,'z',
                   start=c(lon_start,lat_start),
                   count=c(lon_count,lat_count))
# bathy_orig[which(bathy_orig>=0)] <- NA
# bathy_orig <- (-bathy_orig)
hist(bathy_orig)
imagePlot((bathy_orig),asp=1,col=rev(viridis(60)),nlevel=59)
# imagePlot(log10(bathy_orig),asp=1,col=rev(viridis(60)),nlevel=59)

bathy_agg <- aggregate(as.vector(bathy_orig),by=list(lon=lon_c,lat=lat_c),mean,na.rm=T)
bathy_agg <- merge(lonlat_modis,bathy_agg,by=c('lon','lat'),all=T)
bathy_agg_modis <- t(matrix(bathy_agg$x,length(levels(lat_c)),length(levels(lon_c))))
imagePlot(log10(bathy_agg_modis),asp=1,col=rev(viridis(60)),nlevel=59)
imagePlot(log10(bathy_agg_modis),asp=1,col=rev(viridis(60)),nlevel=59)

hist((bathy_agg_modis))
range(bathy_agg_modis,na.rm=T)
min(bathy_agg_modis,na.rm=T)
max(bathy_agg_modis,na.rm=T)
hist((bathy_orig))
range(bathy_orig,na.rm=T)
min(bathy_orig,na.rm=T)
max(bathy_orig,na.rm=T)

setwd('~/Documents/nasa/data/lowres_4km')
# saveRDS(bathy_agg_modis,'gom_wfs_bathy_lowres.rds') # netcdf is smaller and contains metadata
# saveRDS(bathy_orig,'gom_wfs_bathy_hires.rds') # netcdf is smaller and contains metadata


save(bathy_orig,bathy_agg_modis,
     lonlat_bathy,lon_bathy,lat_bathy,
     lonlat_modis,lon_modis,lat_modis,
     file='bathy.RData')
