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

### 3) download sst

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


### SST
url <- 'https://opendap.jpl.nasa.gov/opendap/OceanTemperature/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/2021/001/20210101090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'
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
lonlat_t <- expand.grid(lon=lon_tmp,lat=lat_tmp)
lon_c <- cut(lonlat_t$lon,vec_brk(lon2))
lat_c <- cut(lonlat_t$lat,vec_brk(lat2))
lonlat <- expand.grid(lon=levels(lon_c),lat=levels(lat_c))

times1 <- rep(NA,length(2003:2021))
for(yr in 2003:2021){
  print(paste('Processing',yr, '...',Sys.time()))
  write(paste(Sys.time(), 'Processing MUR SST',yr),'output.txt',append=T)
  ### reference date and julian days
  # yr <- 2021 # 2003:2021
  dates <- data.frame(date=ymd(seq(as.Date(paste0(yr,'-01-01')),as.Date(paste0(yr,'-12-31')),'day')),
                      yday=yday(ymd(seq(as.Date(paste0(yr,'-01-01')),as.Date(paste0(yr,'-12-31')),'day'))))
  
  data_yday <- array(NA,c(length(lon2),
                          length(lat2),
                          nrow(dates)))
  pb <- txtProgressBar(min = 0, max = nrow(dates), style = 3)
  t1 <- system.time(
    for(i in 1:nrow(dates)){
      # url <- 'https://opendap.jpl.nasa.gov/opendap/OceanTemperature/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/2021/001/20210101090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'
      url <- paste0('https://opendap.jpl.nasa.gov/opendap/OceanTemperature/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/',
                    yr,
                    '/',
                    sprintf("%03d",dates$yday[i]),
                    '/',
                    yr,
                    sprintf("%02d",month(dates$date[i])),
                    sprintf("%02d",day(dates$date[i])),
                    '090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc')  
      data <- try(nc_open(url,readunlim=F,suppress_dimvals=T,return_on_error=F))
      if(class(data)!='try-error'){
        time <- ymd_hms(ncatt_get(data,0)$start_time)
        # time <- ncvar_get(data,'time')
        # as.Date(time/(3600*24),origin=as.Date(substr(ncatt_get(data,'time')$units,15,40)))
        
        sst <- ncvar_get(data,'analysed_sst',
                         start=c(lon_start,lat_start,1),
                         count=c(lon_count,lat_count,-1))
        sst <- NISTkTOdegC(sst)
        sst_agg <- aggregate(as.vector(sst),by=list(lon=lon_c,lat=lat_c),mean,na.rm=T)
        sst_agg <- merge(lonlat,sst_agg,by=c('lon','lat'),all=T)
        sst_agg_m <- t(matrix(sst_agg$x,length(levels(lat_c)),length(levels(lon_c))))
        
        nc_close(data)
        data_yday[,,i] <- sst_agg_m
      } else {
        write(paste0(Sys.time(), ' Error MUR SST (i = ',i,') ', url),'output.txt',append=T)
      }
      rm(data,url,sst,sst_agg,sst_agg_m)
      setTxtProgressBar(pb, i)
    }
  )
  write(paste(Sys.time(), 'Processed MUR SST', yr, 'total time (sec):',t1[3]),'output.txt',append=T)
  times1[yr-2002] <- t1[3]
  
  saveRDS(data_yday,paste0('mursst_daily_',yr,'.rds')) # netcdf is smaller and contains metadata
}
cat('total time:',sum(times1,na.rm=T), 'sec')
write(paste('Total time:',sum(times1,na.rm=T), 'sec (2002-2021)'),'output.txt',append=T)
