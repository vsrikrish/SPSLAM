########################################################################################
# compute_rl.R: script to generate ensembles of futures				       #
# 		This file creates an ensemble of maximum local tide gauge	       #
#		observations for each combination of Sea Level Rise and Storm Surge    #
#	        models.								       #
#		It requires the plyr (for manipulating data frames), ncdf4 (for	       #
#		file I/O), and extRemes (for sampling) packages.		       #
#		It also requires BRICK (http://github.com/scrim-network/BRICK).	       #
########################################################################################

# set user controlled parameters
nSOW <- 1e7

year.save <- 2070

burn.in <- 50000	# could be automate dusing G-R diagnostic, but this was chosen heuristically



# load libraries
library(plyr)
library(ncdf4)
library(extRemes)

years <- seq(2020, 2070, by=5)

# find the index of the target year
year.ind <- years - 1850 + 1

# set directories
plot.dir <- 'figures'
data.dir <- 'data'
out.dir <- 'output'

slr.models <- c('linear', 'no accelerated melting', 'accelerated melting', 'accelerated melting')
surge.models <- c('historical', 'stationary', 'stationary', 'nonstationary')

# read in linear fit MCMC
lin.mcmc <- readRDS(file.path(out.dir, 'MCMC-norfolk_linear_slr.rds'))
n.iter <- nrow(lin.mcmc[[1]]$samples)
lin.samp <- lin.mcmc[[1]]$samples[(burn.in+sample(1:nrow(lin.mcmc[[1]]$samples[(burn.in+1):n.iter,]), nSOW, replace=T)),]
colnames(lin.samp) <- c('beta0', 'beta1', 'sigma')

# load and re-fingerprint BRICK ensemble
lat <- 36.96
lon <- -76.33
ncdata <- nc_open(Sys.glob(file.path(data.dir, 'BRICK_physical*.nc')))
slr.gsic <- ncvar_get(ncdata, 'GSIC_RCP85')
slr.gis <- ncvar_get(ncdata, 'GIS_RCP85')
slr.ais <- ncvar_get(ncdata, 'AIS_RCP85')
slr.te <- ncvar_get(ncdata, 'TE_RCP85')
slr.lws <- ncvar_get(ncdata, 'LWS_RCP85')
t.proj <- ncvar_get(ncdata, 'time_proj')
n.ens <- length(ncvar_get(ncdata, 'ens'))
slr.fd <- ncvar_get(ncdata, 'LocalSeaLevel_RCP85') - ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP85')
nc_close(ncdata)
curwd <- setwd('../BRICK/R/')
source('BRICK_LSL.R')
lsl.proj <- brick_lsl(lat.in=lat, lon.in=lon, n.time=length(t.proj), slr_gis=slr.gis, slr_gsic=slr.gsic, slr_ais=slr.ais, slr_te=slr.te, slr_lws=slr.lws)
setwd(curwd)

# load GIA/subsidence ensemble
gia <- read.csv(file.path(data.dir, 'LSLProj_bkgd_299_rcp26.csv'), header=T, colClasses = 'numeric', check.names = F)

# set up 100-year return level storage
rl.100 <- vector('list', length(years))
names(rl.100) <- years

for (i in 1:length(years)) {
  
  SOWs <- vector('list', 4)
  
  # load data from norfolk
  dat <- readRDS(file.path(out.dir, 'processed_norfolk.rds'))
  
  # load BMA weights for mixed surge model
#  weights <- readRDS(file.path(out.dir, 'bma_weights.rds'))
  
  # define list to hold surge model parameters
  surge <- vector('list', 2)
  names(surge) <- c('stationary', 'nonstationary')
  
  # load and postprocess stationary surge model parameters
  # burn in
  burn.in <- 50000
    
  surge.st <- readRDS(file.path(out.dir, 'norfolk_MCMC-stationary.rds'))
  n.iter <- nrow(surge.st[[1]]$samples)
  surge[['stationary']] <- surge.st[[1]]$samples[(burn.in+sample(1:nrow(surge.st[[1]]$samples[(burn.in+1):n.iter,]), nSOW, replace=T)),]
  # load and postprocess stationary surge model parameters
  surge.ns <- readRDS(file.path(out.dir, 'norfolk_MCMC-nonstationary.rds'))
  surge[['nonstationary']] <- surge.ns[[1]]$samples[(burn.in+sample(1:nrow(surge.ns[[1]]$samples[(burn.in+1):n.iter,]), nSOW, replace=T)),]
  
  #######################################################################################################################################
  # sample/calculate SOWs for each case
  # linear/historical record model
  # extrapolate mean sea level
  lin.slr <- lin.samp[,'beta0'] + lin.samp[,'beta1']*(years[i]-dat[['lin']]$year[1]) - dat[['lin']]$mean[which(dat[['lin']]$year == 2015)]
  # bootstrap from the historical storm surge levels
  lin.surge <- sample(dat[['gev']]$ann.max, nSOW, replace=T)
  # combine SLR and surge and convert to m
  SOWs[[1]] <- (lin.slr + lin.surge)/1000
  
  # no fast dynamics SL/stationary surge
  data.brick <- data.frame(LSL=lsl.proj[year.ind[i],] - lsl.proj[(2015-1850+1), ], FD=(slr.fd[year.ind[i],] > 0))
  
  nofd.data <- data.brick[data.brick['FD'] == FALSE,]
  nofd.slr <- sample(nofd.data[,'LSL'], nSOW, replace=T)*1000
  nofd.surge <- revd(nSOW, loc=surge[['stationary']][,'lambda'], 
                     scale=surge[['stationary']][,'sigma'], 
                     shape=surge[['stationary']][,'xi'], 
                     type='GEV')
  if (years[i] %% 10 == 0) {
    nofd.gia <- gia[sample(1:nrow(gia), nSOW, replace=T), as.character(years[i])]*10
  } else {
    nofd.gia <- rowMeans(gia[sample(1:nrow(gia), nSOW, replace=T), as.character(c(floor(years[i]/10)*10, ceiling(years[i]/10)*10))])*10
  }
    
  # combine and convert to m
  SOWs[[2]] <- (nofd.slr + nofd.surge + nofd.gia)/1000
  
  # fast dynamics/stationary surge
  if (sum(data.brick[,'FD']) == 0) {
    fd.data <- nofd.data
  } else {
    fd.data <- data.brick[data.brick['FD'] == TRUE,]
  }
  fdst.slr <- sample(fd.data[,'LSL'], nSOW, replace=T)*1000
  fdst.surge <- revd(nSOW, loc=surge[['stationary']][,'lambda'], 
                     scale=surge[['stationary']][,'sigma'], 
                     shape=surge[['stationary']][,'xi'], 
                     type='GEV')
  if (years[i] %% 10 == 0) {
    fdst.gia <- gia[sample(1:nrow(gia), nSOW, replace=T), as.character(years[i])]*10
  } else {
    fdst.gia <- rowMeans(gia[sample(1:nrow(gia), nSOW, replace=T), as.character(c(floor(years[i]/10)*10, ceiling(years[i]/10)*10))])*10
  }
  SOWs[[3]] <- (fdst.slr + fdst.surge + fdst.gia)/1000
  
  # fast dynamics/BMA surge
  fdbma.slr <- sample(fd.data[,'LSL'], nSOW, replace=T)*1000
  
#  bma.case <- sample(1:2, nSOW, replace=T, prob=weights)
  temp <- dat[['temps']]$forcing[dat[['temps']]$forcing$year == years[i],'temp']
#  idx <- sample(1:nrow(surge[['stationary']]), nSOW, replace=T)
  # lambda <- numeric(nSOW)
  # sigma <- numeric(nSOW)
  # xi <- numeric(nSOW)
  # 
  # lambda[bma.case == 1] <- surge[['stationary']][idx[bma.case == 1], 'lambda']
  # lambda[bma.case == 2] <- surge[['nonstationary']][idx[bma.case == 2], 'lambda0'] + surge[['nonstationary']][idx[bma.case == 2], 'lambda1'] * temp.2065
  # sigma[bma.case == 1] <- surge[['stationary']][idx[bma.case == 1], 'sigma']
  # sigma[bma.case == 2] <- surge[['nonstationary']][idx[bma.case == 2], 'sigma']
  # xi[bma.case == 1] <- surge[['stationary']][idx[bma.case == 1], 'xi']
  # xi[bma.case == 2] <- surge[['nonstationary']][idx[bma.case == 2], 'xi']
  
#  fdbma.surge <- revd(nSOW, loc=lambda, scale=sigma, shape=xi, type='GEV')
  fdbma.surge <- revd(nSOW, 
                      loc=surge[['nonstationary']][,'lambda0'] + surge[['nonstationary']][,'lambda1'] * temp,
                      scale=surge[['nonstationary']][,'sigma'],
                      shape=surge[['nonstationary']][,'xi'],
                      type='GEV')
  if (years[i] %% 10 == 0) {
    fdbma.gia <- gia[sample(1:nrow(gia), nSOW, replace=T), as.character(years[i])]*10
  } else {
    fdbma.gia <- rowMeans(gia[sample(1:nrow(gia), nSOW, replace=T), as.character(c(floor(years[i]/10)*10, ceiling(years[i]/10)*10))])*10
  }
  SOWs[[4]] <- (fdbma.slr + fdbma.surge + fdbma.gia)/1000
  #######################################################################################################################
  # save SOWsbma.case <- sample(1:2, nSOW, replace=T, prob=weights)
  
  if (years[i] == year.save) {
    saveRDS(SOWs, file.path(out.dir, 'SL-SOWs.rds'))
  }

  # compute 100-year return levels
  rl.100[[i]] <- lapply(SOWs, quantile, prob=c(.99))
}

saveRDS(rl.100, file.path(out.dir, 'rl100.rds'))
