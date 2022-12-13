### modis derived products and netcdf creation

library(abind)
library(fields)
library(lubridate)
library(ncdf4)
library(ncdf4.helpers)
library(NISTunits)

vec_brk <- function (input) {
  core <- (diff(input)/2) + input[1:(length(input)-1)]
  beg <- input[1] - (diff(input[1:2])/2)
  last <- input[length(input)] + (diff(input[(length(input)-1):(length(input))])/2)
  output <- c(beg,core,last)
  output
}

rrs_nlw <- function (data, wavelength, parm_vec = parm, F0_vec = sat_wavelength){
  pinx <- grep(wavelength, parm_vec)
  finx <- which(F0_vec==wavelength)
  return(data[pinx,,,] * F0[finx])
}

### 1) download daily 4km Aqua data
# a) what bounding box
# (-87.5 30.7, -81, 24.2)
lonbox_w <- -87.5 ### mouth of Mississippi River
latbox_n <- 30.7 ### northern coast
lonbox_e <- -81 ### Florida Bay
latbox_s <- 24.2 ### southern edge of Key West

# Rrs = nLw/F0
# Rrs: sr^-1; nLw: mW cm^-2 um^-1 sr^-1; F0: mW cm^-2 um^-1
sat_wavelength <- c(412,443,469,488,531,547,551,555,645,667,678,748,859,869,1240,1640,2130)
F0 <- c(172.912,187.622,205.878,194.933,185.747,186.539,186.539,183.869,157.811,152.255,148.052,128.065,97.174,95.824,45.467,23.977,9.885)

parms <- c('CHL_chlor_a','FLH_nflh','RRS_Rrs_443','RRS_Rrs_488','RRS_Rrs_531','RRS_Rrs_547','RRS_Rrs_555','RRS_Rrs_667','RRS_Rrs_678')
parm <- substr(parms,5,11)

# url <- 'http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2011/002/A2011002.L3m_DAY_RRS_Rrs_443_4km.nc'
url <- 'http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2011/0102/AQUA_MODIS.20110102.L3m.DAY.RRS.Rrs_443.4km.nc'
system.time(modis1 <- nc_open(url,readunlim=F,suppress_dimvals=T,return_on_error=T))
# atts <- ncatt_get(modis1,0)
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
nc_close(modis1)
rm(modis1)


setwd('~/Documents/nasa/data/lowres_4km')
# yr <- 2003
pb <- txtProgressBar(min = 0, max = length(2003:2021), style = 3)
t1 <- system.time(
  for(yr in 2003:2021){
    # data_yday <- readRDS(paste0('modisa_daily_',yr,'.rds'))
    data_yday <- readRDS(paste0('aqua_modis_daily_',yr,'_4km.rds'))
    data_out <- array(NA,c(12,
                           length(lon2),
                           length(lat2),
                           dim(data_yday)[4]))
    if(yr>2002){
      data_previous <- readRDS(paste0('aqua_modis_daily_',(yr-1),'_4km.rds'))
      begin <- dim(data_previous)[4]-74
      last <- dim(data_previous)[4]
      data_previous <- data_previous[1:2,,,begin:last]
    }
    
    # nlw_443
    nlw_443 <- rrs_nlw(data_yday, 443)
    
    # nlw_488
    nlw_488 <- rrs_nlw(data_yday, 488)
    
    # nlw_531
    nlw_531 <- rrs_nlw(data_yday, 531)
    
    # nlw_667
    nlw_667 <- rrs_nlw(data_yday, 667)
    
    # nlw_678
    nlw_678 <- rrs_nlw(data_yday, 678)
    
    # ssnlw488 - nlw_443, nlw_488, nlw_531
    ssnlw488 <- nlw_488 - nlw_443 - (nlw_531 - nlw_443) * (488 - 443) / (531 - 443)
    
    # RBD - nlw_667, nlw_678
    rbd <- nlw_678 - nlw_667
    
    # KBBI - nlw_667, nlw_678
    kbbi <- rbd / (nlw_678 + nlw_667)
    
    # ABI - nflh, rrs_547
    # nflh / (1 + (rrs_547 - 0.0015) * 80)
    abi <- data_yday[grep('nflh', parm),,,] / (1 + (data_yday[grep(547, parm),,,] - 0.0015) * 80)
    
    # carder_bbp - rrs_555
    # (2.058 * Rrs_555 - .00182)
    carder_bbp <- 2.058 * data_yday[grep('555', parm),,,] - .00182
    
    # morel_bbp - chlor_a
    # (.3 * pow(chlor_a, .62) * (.002 + .02 * (.5 - .25 * log10(chlor_a))))
    morel_bbp <- .3 * (data_yday[grep('chlor_a', parm),,,]^.62) * (.002 + .02 * (.5 - .25 * log10(data_yday[grep('chlor_a', parm),,,])))
    
    ### Carder to Morel ratio
    cm_bbp <- carder_bbp / morel_bbp
    
    # chlor_a anomaly
    chl_yday <- data_yday[grep('chlor_a', parm),,,]
    chl_yday <- abind(data_previous[grep('chlor_a', parm),,,],chl_yday,along=3)
    chl_anom2 <- array(NA,c(length(lon2),length(lat2),dim(data_yday)[4]))
    # nflh anomaly
    nflh_yday <- data_yday[grep('nflh', parm),,,]
    nflh_yday <- abind(data_previous[grep('nflh', parm),,,],nflh_yday,along=3)
    nflh_anom2 <- array(NA,c(length(lon2),length(lat2),dim(data_yday)[4]))
    n <- 1 #counter
    for(k in 76:dim(chl_yday)[3]){
      ### 15 day lagged anomaly of 60 day mean
      # chlor_A
      chl_lm_60d <- apply(chl_yday[,,(k-75):(k-16)],c(1,2),mean,na.rm=T)
      chl_anom2[,,n] <- chl_yday[,,k] - chl_lm_60d
      # nflh
      nflh_lm_60d <- apply(nflh_yday[,,(k-75):(k-16)],c(1,2),mean,na.rm=T)
      nflh_anom2[,,n] <- nflh_yday[,,k] - nflh_lm_60d
      n <- n + 1
      rm(chl_lm_60d,nflh_lm_60d)
    }
    
    data_out[1,,,] <- data_yday[grep('chlor_a', parm),,,]
    data_out[2,,,] <- chl_anom2
    data_out[3,,,] <- data_yday[grep('nflh', parm),,,]
    data_out[4,,,] <- nflh_anom2
    data_out[5,,,] <- data_yday[grep('Rrs_667', parm),,,]
    data_out[6,,,] <- ssnlw488
    data_out[7,,,] <- carder_bbp
    data_out[8,,,] <- morel_bbp
    data_out[9,,,] <- cm_bbp
    data_out[10,,,] <- abi
    data_out[11,,,] <- rbd
    data_out[12,,,] <- kbbi
    
    saveRDS(data_out,paste0('aqua_modis_daily_input_',yr,'_4km.rds')) # netcdf is smaller and contains metadata 
    rm(data_yday,data_previous,begin,last,data_out,chl_anom2,chl_yday,nflh_anom2,nflh_yday,abi,carder_bbp,morel_bbp,cm_bbp,ssnlw488,rbd,kbbi,nlw_443,nlw_488,nlw_531,nlw_667,nlw_678)
    # gc()
    setTxtProgressBar(pb, yr-2002)
  }
)
t1

### what are the output variables?
parm_out <- c('chlor_a','chl_anom','nflh','nflh_anom','rrs_667','ssnlw488','carder_bbp','morel_bbp','cm_bbp','abi','rbd','kbbi')


### updated code 2022/12/12 to fix nflh_anom calculation and include cm_bbp
# ### Carder to Morel ratio
# 
# pb <- txtProgressBar(min = 0, max = length(2003:2021), style = 3)
# t1 <- system.time(
#   for(yr in 2003:2021){
#     setwd('~/Documents/nasa/data/lowres_4km')
#     # data_yday <- readRDS(paste0('modisa_daily_',yr,'.rds'))
#     data_yday <- readRDS(paste0('aqua_modis_daily_input_',yr,'_4km.rds'))
#     
#     cm_bbp <- data_yday[grep('carder',parm_out),,,]/data_yday[grep('morel',parm_out),,,]
#     
#     new_data <- abind(data_yday,cm_bbp,along=1)
#     setwd('~/Documents/nasa/data')
#     saveRDS(new_data,paste0('aqua_modis_daily_input_',yr,'_4km.rds')) # netcdf is smaller and contains metadata
#     
#     setTxtProgressBar(pb, yr-2002)
#     rm(data_yday,cm_bbp,new_data)
#   }
# )
# 
# parm_out <- c('chlor_a','chl_anom','nflh','nflh_anom','rrs_667','abi','bbp_carder','bbp_morel','ssnlw488','rbd','kbbi','cm_bbp')
