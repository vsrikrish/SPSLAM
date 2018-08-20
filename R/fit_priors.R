##################################################################################
# fit_priors.R: This file fits prior distributions based on ML estimates	 #
# 		at other tide gauge stations.					 #
#		It requires the R packages DEoptim (for finding MLEs) and	 #
#		ncdf4 (to store output).					 #
##################################################################################

library(DEoptim)
library(ncdf4)
source('~/work/ScienceUncertainty/R/gev_posterior.R')

fit_priors_stations <- function(tide.data.path, temp.data.path, model) {
  
  print('starting to process many tide gauge stations data')
  
  # read in tide gauge data
  
  filetype='csv'
  files.tg <- list.files(path=tide.data.path, pattern=filetype)
  
  data_set <- vector('list', length(files.tg))
  for (dd in 1:length(files.tg)) {
    # print an update of progress to the screen
    print(paste('now reading in data set ',dd,' / ',length(files.tg),sep=''))
    names(data_set)[dd] <- substr(files.tg[dd], start=1, stop=7)
    data.tmp <- read.table(file.path(tide.data.path,files.tg[dd]), header=FALSE, sep=',')
    data_set[[dd]] <- data.frame(matrix(ncol=2, nrow=nrow(data.tmp)))
    colnames(data_set[[dd]]) <- c('year','sl')
    data_set[[dd]]$year <- data.tmp$V1
    data_set[[dd]]$sl <- data.tmp$V5
  }
  
  # compute annual SL means
  sl.means <- lapply(data_set, function(l) aggregate(l$sl, list(l$year), mean))
  sl.means <- lapply(sl.means, setNames, c('year', 'mean'))
  
  # compute residuals from means
  dat.combined <- lapply(1:length(sl.means), function(i) merge(data_set[[i]], sl.means[[i]]))
  dat.combined <- lapply(dat.combined, function(l) cbind(l, l$sl - l$mean))
  dat.combined <- lapply(dat.combined, setNames, c('year', 'sl', 'mean', 'resid'))
  
  # compute annual block maxima
  data_many <- lapply(dat.combined, function(l) aggregate(l$resid, list(l$year), max))
  names(data_many) <- names(data_set)
  data_many <- lapply(data_many, setNames, c('year', 'ann.max'))
  
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
  
  temperature_df <- data.frame(year=time_forc, temp=temperature_forc)
  
  if (model == 'stationary') {
    parnames <- c('lambda', 'sigma', 'xi')
  } else if (model == 'nonstationary') {
    parnames <- c('lambda0', 'lambda1', 'sigma', 'xi')
  }
  
  # set DEOptim parameters
  NP.deoptim <- 100
  # number of DEoptim iterations
  niter.deoptim <- 200
  # differential weighting factor
  F.deoptim <- 0.8
  # crossover probability
  CR.deoptim <- 0.9
  # lower and upper bounds for each variable
  if (model == 'stationary') {
    bounds.lower <- c(0, 0, -2)
    bounds.upper <- c(5000, 1000, 2)
  } else if (model == 'nonstationary') {
    bounds.lower <- c(0, 0, 0, -2)
    bounds.upper <- c(5000, 2000, 1000, 2)
  }
  
  # compute MLE estimates for each site
  est.ml <- data.frame(matrix(ncol=length(parnames), nrow=length(data_many)))
  colnames(est.ml) <- parnames
  for (dd in 1:length(data_many)) {
    print(paste('finding MLE for', dd, '/', length(data_many), sep=' '))
    # merge SLR with temperatures
    dat <- merge(data_many[[dd]], temperature_df)
    colnames(dat) <- c('year', 'ann.max', 'temp')
    
    if (model == 'stationary') {
      aux_data <- NULL
    } else if (model == 'nonstationary') {
      aux_data <- dat$temp
    }
    deoptim_out <- DEoptim(fn=neg_log_likelihood, 
                           lower=bounds.lower,
                           upper=bounds.upper,
                           DEoptim.control(NP=NP.deoptim, itermax=niter.deoptim, F=F.deoptim, CR=CR.deoptim, trace=FALSE),
                           parnames=parnames, model=model, dat=dat$ann.max, aux_data=aux_data)

    est.ml[dd,] <- deoptim_out$optim$bestmem
  }
  
  # fit distributions: normal for any location and shape parameters, gamma for scale (to ensure positivity)
  dist.means <- colMeans(est.ml)
  dist.sd <- apply(est.ml, 2, sd)
  
  # save fit moments as list
  priors <- vector('list', length(parnames))
  names(priors) <- parnames
  for (par in parnames) {
    priors[[par]] <- vector('list', 3)
    if (par == 'sigma') {
      names(priors[[par]]) <- c('type', 'shape', 'scale')
      priors[[par]]$type <- 'gamma'
      priors[[par]]$scale <- dist.sd[match('sigma', parnames)]^2/dist.means[match('sigma', parnames)]
      priors[[par]]$shape <- dist.means[match('sigma', parnames)]/priors[[par]]$scale
    } else {
      names(priors[[par]]) <- c('type', 'mean', 'sd')
      priors[[par]]$type <- 'normal'
      priors[[par]]$mean <- dist.means[match(par, parnames)]
      priors[[par]]$sd <- dist.sd[match(par, parnames)]
    }
  }
  
  # return prior list
  priors
  
}