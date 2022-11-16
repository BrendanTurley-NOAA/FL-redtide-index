### modis derived products and netcdf creation

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


## red band difference
### Ruhul Amin, Jing Zhou, Alex Gilerson, Barry Gross, Fred Moshary, and Samir Ahmed, "Novel optical techniques for detecting and classifying toxic dinoflagellate Karenia brevis blooms using satellite imagery," Opt. Express 17, 9126-9144 (2009)
# nLw678 - nLw667
# (nLw678 - nLw667)/(nLw678 + nLw667)

### Karenia brevis bloom index (KBBI)


# Rrs = nLw/F0
# Rrs: sr^-1; nLw: mW cm^-2 um^-1 sr^-1; F0: mW cm^-2 um^-1
sat_wavelength <- c(412,443,469,488,531,547,551,555,645,667,678,748,859,869,1240,1640,2130)
F0 <- c(172.912,187.622,205.878,194.933,185.747,186.539,186.539,183.869,157.811,152.255,148.052,128.065,97.174,95.824,45.467,23.977,9.885)

parms <- c('CHL_chlor_a','FLH_nflh','RRS_Rrs_443','RRS_Rrs_488','RRS_Rrs_531','RRS_Rrs_547','RRS_Rrs_555','RRS_Rrs_667','RRS_Rrs_678')
parm <- substr(parms,5,11)

url <- 'http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2011/002/A2011002.L3m_DAY_RRS_Rrs_443_4km.nc'
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


setwd('~/Documents/nasa/data/lowres_4km')
yr <- 2003
data_yday <- readRDS(paste0('modisa_daily_',yr,'.rds'))

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

# bbp_Morel - chlor_a
# (.3 * pow(chlor_a, .62) * (.002 + .02 * (.5 - .25 * log10(chlor_a))))
bbp_morel <- .3 * (data_yday[grep('chlor_a', parm),,,]^.62) * (.002 + .02 * (.5 - .25 * log10(data_yday[grep('chlor_a', parm),,,])))

# bbp_Carder - rrs_555
# (2.058 * Rrs_555 - .00182)
bbp_Carder <- 2.058 * data_yday[grep('555', parm),,,] - .00182

# chlor_a anomaly
chl_yday <- data_yday[grep('chlor_a', parm),,,]
chl_anom <- array(NA,c(length(lon2),length(lat2),length(75:365)))
for(i in 1:length(lon2)){
  for(j in 1:length(lat2)){
    n <- 1
    for(k in 75:365){
      ### 15 day lagged anomaly of 60 day mean
      lm_60d <- mean(chl_yday[i,j,(k-74):(k-74+59)],na.rm=T)
      chl_anom[i,j,n] <- chl_yday[i,j,k] - lm_60d
      n <- n + 1
    }
  }
}

dim(chl_yday[1,,])


# chlor_a anomaly
chl_yday <- data_yday[grep('chlor_a', parm),,,]
chl_anom2 <- array(NA,c(length(lon2),length(lat2),length(75:365)))
n <- 1
for(k in 75:365){
  ### 15 day lagged anomaly of 60 day mean
  lm_60d <- apply(chl_yday[,,(k-74):(k-74+59)],c(1,2),mean,na.rm=T)
  chl_anom2[,,n] <- chl_yday[,,k] - lm_60d
  n <- n + 1
}

identical(chl_anom,chl_anom2)


