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

### 1) download daily 4km Aqua data
# a) what bounding box
# (-87.5 30.7, -81, 24.2)
lonbox_w <- -87.5 ### mouth of Mississippi River
latbox_n <- 30.7 ### northern coast
lonbox_e <- -81 ### Florida Bay
latbox_s <- 24.2 ### southern edge of Key West

# b) which parameters
# chlor_a
# nflh
# rrs_443
# rrs_488
# rrs_531
# rrs_547
# rrs_555
# rrs_667
# rrs_678

# c) what derived parameters
# ABI - rrs_547 & nflh
# bbp_Morel - chlor_a
# bbp_Carder - rrs_555
# ssnlw488 - nlw_443, nlw_488, nlw_531
# RBD - nlw_667, nlw_678
# KBBI - nlw_667, nlw_678

## red band difference
### Ruhul Amin, Jing Zhou, Alex Gilerson, Barry Gross, Fred Moshary, and Samir Ahmed, "Novel optical techniques for detecting and classifying toxic dinoflagellate Karenia brevis blooms using satellite imagery," Opt. Express 17, 9126-9144 (2009)
# nLw678 - nLw667
### Karenia brevis bloom index (KBBI)
# (nLw678 - nLw667)/(nLw678 + nLw667)
# Rrs = nLw/F0
# Rrs: sr^-1; nLw: mW cm^-2 um^-1 sr^-1; F0: mW cm^-2 um^-1
sat_wavelength <- c(412,443,469,488,531,547,551,555,645,667,678,748,859,869,1240,1640,2130)
F0 <- c(172.912,187.622,205.878,194.933,185.747,186.539,186.539,183.869,157.811,152.255,148.052,128.065,97.174,95.824,45.467,23.977,9.885)

parms <- c('CHL.chlor_a','FLH.nflh','RRS.Rrs_443','RRS.Rrs_488','RRS.Rrs_531','RRS.Rrs_547','RRS.Rrs_555','RRS.Rrs_667','RRS.Rrs_678')
parm <- substr(parms,5,11)

# url <- 'http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2011/002/A2011002.L3m_DAY_RRS_Rrs_443_4km.nc'
url <- 'http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2011/0102/AQUA_MODIS.20110102.L3m.DAY.RRS.Rrs_443.4km.nc'
# system.time(modis1 <- nc_open(url,readunlim=T,suppress_dimvals=F,return_on_error=T))
system.time(modis1 <- nc_open(url,readunlim=F,suppress_dimvals=T,return_on_error=T))
atts <- ncatt_get(modis1,0)
# names(atts)[-c(2:8,10:11,16:18,27:28,31:35,38,43:47,49:50,52,54,56:60)]
# names(atts)[c(2:8,10:11,16:18,27:28,31:35,38,43:47,49:50,52,54,56:60)]
att_ind <- seq(1,length(atts))[-c(2:8,10:11,16:18,27:28,31:35,38,43:47,49:50,52,54,56:60)]
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

times1 <- rep(NA,length(2002:2021))
for(yr in 2002:2021){ # 2022-11-08; 2003-2016,2021 completed
  print(paste('Processing',yr, '...',Sys.time()))
  write(paste(Sys.time(), 'Processing',yr),'output.txt',append=T)
  ### reference date and julian days
  # yr <- 2021 # 2002:2021
  mth <- ifelse(yr==2002,11,01)
  dd <- ifelse(yr==2002,15,01)
  dates <- data.frame(date=ymd(seq(as.Date(paste0(yr,'-',sprintf("%02d",mth),'-',sprintf("%02d",dd))),as.Date(paste0(yr,'-12-31')),'day')),
                      yday=yday(ymd(seq(as.Date(paste0(yr,'-',sprintf("%02d",mth),'-',sprintf("%02d",dd))),as.Date(paste0(yr,'-12-31')),'day'))))
  ## create netcdf file
  dimlon <- ncdim_def('Lon','degreesE',lon2)
  dimlat <- ncdim_def('Lat','degreesN',lat2)
  dates1 <- as.POSIXct(paste0(dates$date,00:00),tz='GMT')
  dates2 <- as.numeric(dates1)/86400
  # as.Date(dates2[1],origin='1970-01-01')
  dimtime <- ncdim_def('Time','days since 1970-01-01',dates2)
  chlor_a <- ncvar_def('chlor_a','mg m^-3',list(dimlon,dimlat,dimtime),-32767,"Chlorophyll Concentration, OCI Algorithm",prec='double',compression=5) # prec='float may be smaller'
  nflh <- ncvar_def('nflh','W m^-2 um^-1 sr^-1',list(dimlon,dimlat,dimtime),-32767,"Normalized Fluorescence Line Height",prec='double',compression=5) # prec='float may be smaller'
  rrs_443 <- ncvar_def('Rrs_443','sr^-1',list(dimlon,dimlat,dimtime),-32767,"Remote sensing reflectance at 443 nm",prec='double',compression=5) # prec='float may be smaller'
  rrs_488 <- ncvar_def('Rrs_488','sr^-1',list(dimlon,dimlat,dimtime),-32767,"Remote sensing reflectance at 488 nm",prec='double',compression=5) # prec='float may be smaller'
  rrs_531 <- ncvar_def('Rrs_531','sr^-1',list(dimlon,dimlat,dimtime),-32767,"Remote sensing reflectance at 531 nm",prec='double',compression=5) # prec='float may be smaller'
  rrs_547 <- ncvar_def('Rrs_547','sr^-1',list(dimlon,dimlat,dimtime),-32767,"Remote sensing reflectance at 547 nm",prec='double',compression=5) # prec='float may be smaller'
  rrs_555 <- ncvar_def('Rrs_555','sr^-1',list(dimlon,dimlat,dimtime),-32767,"Remote sensing reflectance at 555 nm",prec='double',compression=5) # prec='float may be smaller'
  rrs_667 <- ncvar_def('Rrs_667','sr^-1',list(dimlon,dimlat,dimtime),-32767,"Remote sensing reflectance at 667 nm",prec='double',compression=5) # prec='float may be smaller'
  rrs_678 <- ncvar_def('Rrs_678','sr^-1',list(dimlon,dimlat,dimtime),-32767,"Remote sensing reflectance at 678 nm",prec='double',compression=5) # prec='float may be smaller'
  # modis_tmp <- nc_create('modis_tmp.nc',list(chlor_a,nflh,rrs_443,rrs_488,rrs_531,rrs_547,rrs_555,rrs_667,rrs_678))
  modis_tmp <- nc_create(paste0('modisa_daily_',yr,'.nc'),list(chlor_a,nflh,rrs_443,rrs_488,rrs_531,rrs_547,rrs_555,rrs_667,rrs_678))
  nc.copy.atts(modis1,0,modis_tmp,0,names(atts)[att_ind]) # not tested
  ncatt_put(modis_tmp,0,"northernmost_latitude",max(lat2),prec='float')
  ncatt_put(modis_tmp,0,"southernmost_latitude",min(lat2),prec='float')
  ncatt_put(modis_tmp,0,"westernmost_longitude",min(lon2),prec='float')
  ncatt_put(modis_tmp,0,"easternmost_longitude",max(lon2),prec='float')
  ncatt_put(modis_tmp,0,"time_coverage_start",paste(dates1[1]))
  ncatt_put(modis_tmp,0,"time_coverage_end",paste(dates1[length(dates1)]))
  ncatt_put(modis_tmp,0,"modified_by",'Brendan Turley')
  
  data_yday <- array(NA,c(9,
                          length(lon2),
                          length(lat2),
                          nrow(dates)))
  pb <- txtProgressBar(min = 0, max = nrow(dates), style = 3)
  t1 <- system.time(
    for(i in 1:nrow(dates)){
      for(j in 1:length(parms)){
        url <- paste0('http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/',
                      yr,
                      '/',
                      paste0(sprintf("%02d",month(dates$date[i])),
                             sprintf("%02d",day(dates$date[i]))),
                      '/AQUA_MODIS.',
                      yr,
                      paste0(sprintf("%02d",month(dates$date[i])),
                             sprintf("%02d",day(dates$date[i]))),
                      '.L3m.DAY.',
                      parms[j],
                      '.4km.nc')
        # url <- paste0('http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/',
        #               yr,
        #               '/',
        #               sprintf("%03d",dates$yday[i]),
        #               '/A',
        #               yr,
        #               sprintf("%03d",dates$yday[i]),
        #               '.L3m_DAY_',
        #               parms[j],
        #               '_4km.nc')
        # modis <- try(nc_open(url))
        modis <- try(nc_open(url,readunlim=F,suppress_dimvals=T,return_on_error=F))
        if(class(modis)!='try-error'){
          data <- ncvar_get(modis,parm[j],start=c(lon_start,lat_start),count=c(lon_count,lat_count))
          # tlon <- ncvar_get(modis, 'lon',start=c(lon_start),count=c(lon_count))
          # tlat <- ncvar_get(modis, 'lat',start=c(lat_start),count=c(lat_count))
          if(i==1){
            nc.copy.atts(modis,parm[j],modis_tmp,parm[j],c('name','units','dim','_FillValue','longname','prec')) # not tested
          }
          nc_close(modis)
          data_yday[j,,,i] <- data
        } else {
          write(paste0(Sys.time(), ' Error (i = ',i,', j = ',j,') ', url),'output.txt',append=T)
        }
        rm(modis,data,url)
      }
      setTxtProgressBar(pb, i)
    }
  )
  write(paste(Sys.time(), 'Processed', yr, 'total time (sec):',t1[3]),'output.txt',append=T)
  times1[yr-2002] <- t1[3]
  # "2022-11-04 11:26:35 CDT"
  # user   system  elapsed 
  # 134.101   48.079 4389.390 
  
  ### fill netcdf file
  ncvar_put(modis_tmp,chlor_a,data_yday[1,,,])
  ncvar_put(modis_tmp,nflh,data_yday[2,,,])
  ncvar_put(modis_tmp,rrs_443,data_yday[3,,,])
  ncvar_put(modis_tmp,rrs_488,data_yday[4,,,])
  ncvar_put(modis_tmp,rrs_531,data_yday[5,,,])
  ncvar_put(modis_tmp,rrs_547,data_yday[6,,,])
  ncvar_put(modis_tmp,rrs_555,data_yday[7,,,])
  ncvar_put(modis_tmp,rrs_667,data_yday[8,,,])
  ncvar_put(modis_tmp,rrs_678,data_yday[9,,,])
  ncatt_put(modis_tmp,0,"last_Modified",paste0(with_tz(Sys.time(),tz='utc'),'Z'))
  nc_close(modis_tmp)
  
  setwd('~/Documents/nasa/data/lowres_4km')
  saveRDS(data_yday,paste0('aqua_modis_daily_',yr,'_4km.rds')) # netcdf is smaller and contains metadata
  # data_yday <- readRDS('modisa_daily_2021.rds')
}
cat('total time:',sum(times1,na.rm=T), 'sec')
# write(paste('Total time:',sum(times1,na.rm=T), 'sec (2017-2020)'),'output.txt',append=T)
