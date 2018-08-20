# script to process data for annual block maxima for Norfolk, VA, USA

# produces the list object 'data_norfolk', which has all of the necessary 
# information to estimate the GEV parameters via MLE
#
# relative paths are ok, but absolute paths are probably better

library(ncdf4)

# function to process data files from Norfolk
data_process_norfolk <- function(sl.data.path, temp.data.path) {
  
  
  # create the object to hold the calibration information for the Norfolk site
  # holds block maximum and global mean temperature information
  data_norfolk <- vector('list', 3)
  names(data_norfolk) <- c('gev', 'lin', 'temps')
  
  
  # read in tide gauge data (hourly series)
  
  filetype <- 'txt'
  septype <- '\t'
  
  files.tg <- list.files(path=sl.data.path,pattern=filetype)
  
  dat <- read.table(file.path(sl.data.path,files.tg[1]), header = TRUE, sep=septype)
  if(length(files.tg) > 1) {
    for (ff in 2:length(files.tg)) {
      dat <- rbind(dat, read.table(file.path(sl.data.path,files.tg[ff]), header = TRUE, sep=septype))
    }
  }
  
  # convert sea levels from m to mm, consistent with the other data
  dat$sl <- 1000* dat$sl
  
  # get year from date string
  dat$year   <- as.numeric(substr(as.character(dat$date), start=1, stop=4))
  
  # compute annual SL means
  sl.means <- aggregate(dat$sl, list(dat$year), mean)
  colnames(sl.means) <- c('year', 'mean')
  data_norfolk[['lin']] <- sl.means
  
  # compute residuals from means
  dat.combined <- merge(dat, sl.means)
  dat.combined$resid <- dat.combined$sl - dat.combined$mean
  
  # compute block maxima
  data_norfolk[['gev']] <- aggregate(dat.combined$resid, list(dat.combined$year), max)
  colnames(data_norfolk[['gev']]) <- c('year', 'ann.max')
  
  print('done processing the Sewells Point/Norfolk tide gauge data set')
  
  ## process global mean temperature data
  print('starting to process global mean temperatures')
  
  # read the projections temperature forcing from CRNM (CMIP5)
  # note: these are in K, but going to normalize, so will take a difference and
  # it's same as celsius
  ncdata <- nc_open(file.path(temp.data.path, 'global.tas.aann.CNRM-CM5.historical+rcp85.r1i1p1.18500101-21001231.nc'))
  temperature_proj <- ncvar_get(ncdata, 'tas')
  time_proj <- ncvar_get(ncdata, 'time')
  nc_close(ncdata)
  
  # 'time' on the netcdf ile is YYYYMMDD, where MMDD is July 2 each year (becuase
  # of the averaging). so pluck off just the years
  time_proj <- floor(time_proj/10000)
  
  # read the historical forcing from NOAA
  data.tmp <- read.table(file.path(temp.data.path, 'noaa_temperature_1880-2017.csv'), header = TRUE, sep=',')
  time_hist <- data.tmp$Year
  temperature_hist <- data.tmp$Value

  data.tmp <- read.table(file.path(temp.data.path, 'HadCRUT.4.4.0.0.annual_ns_avg.txt'))
  time_hadcrut <- data.tmp[,1]
  temperature_hadcrut <- data.tmp[,2]  
  
  ind_norm <- which(time_hist==1901):which(time_hist==2000)
  temperature_hist <- temperature_hist - mean(temperature_hist[ind_norm])
  
  ind_norm <- which(time_hadcrut==1901):which(time_hadcrut==2000)
  temperature_hadcrut <- temperature_hadcrut - mean(temperature_hadcrut[ind_norm])
  
  ind_norm <- which(time_proj==1901):which(time_proj==2000)
  temperature_proj <- temperature_proj - mean(temperature_proj[ind_norm])
  
  # set up the forcing, wherein the historical starts, then projections finish
  # 1850-1880 -- hadcrut4
  # 1880-2016 -- NOAA historical
  # 2017-2100 -- CRNM projection
  time_forc <- min(time_hadcrut):max(time_proj)
  temperature_forc <- rep(NA, length(time_forc))
  
  ind_hadcrut <- which(time_hadcrut==time_forc[1]):(which(time_hadcrut==time_hist[1])-1)
  ind_hist    <- 1:length(time_hist)
  ind_proj    <- which(time_proj==(max(time_hist)+1)):which(time_proj==max(time_forc))
  
  temperature_forc <- c(temperature_hadcrut[ind_hadcrut],
                        temperature_hist[ind_hist]      ,
                        temperature_proj[ind_proj]      )
  
  # maximum temperature serves as an additinoal prior constraint on kappa0, kappa1
  # that is, kappa1 > -kappa0/Tmax (otherwise, kappa = kappa0 + kappa1*T might be
  # negative)
  Tmax <- max(temperature_forc)
  
  data_norfolk[['temps']] <- vector('list', 2)
  names(data_norfolk[['temps']]) <- c('forcing', 'max')
  data_norfolk[['temps']]$forcing <- data.frame(year=time_forc, temp=temperature_forc)
  data_norfolk[['temps']]$max <- Tmax
  
  data_norfolk

}
