##############################################################################################
# calibrate_model.R: main file for MCMC calibration for GEV models			     #
# 		     This file requires the parallel and adaptMCMC packages and will	     #
#		     	  need some reworking to run in serial.				     #
##############################################################################################


library(parallel)
library(adaptMCMC)

# set paths for data and storage
data.path <- '~/work/ScienceUncertainty/data'
output.path <- '~/work/ScienceUncertainty/output'
R.path <- '~/work/ScienceUncertainty/R'

# create directories if they don't already exist
for (path in c(data.path, output.path)) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

args <- commandArgs(trailingOnly = TRUE)

source(file.path(R.path, 'gev_posterior.R'))

# set parameters for this run
# method can be "MLE," "prelim," or "production"
method <- args[1]
# model can be 'stationary' or 'nonstationary' (location parameter is a linear function of temp)
model <- args[2]

# these are required for either MCMC method
if (method == 'production') {
  # number of MCMC iterations per chain
  n.iter <- 5e5
  # number of MCMC chains
  n.chain <- 10
  # number of cpus
  n.cpu <- parallel::detectCores()
} else if (method == 'prelim') {
  n.iter <- 1e4
  n.chain <- 1
  n.cpu <- 1
} else {
# these are required for MLE
  library(DEoptim)
  # DEOptim control parameters
  # number of DEOptim population members
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
    bounds.upper <- c(2000, 800, 2)
  } else if (model == 'nonstationary') {
    bounds.lower <- c(0, 0, 0, -2)
    bounds.upper <- c(2000, 500, 800, 2)
  }
}

# process data
filename.data <- 'processed_norfolk.rds'
# don't need to re-process if this has already been done
if (!file.exists(file.path(output.path, filename.data))) {
  source(file.path(R.path, 'process_norfolk_data.R'))
  dat_list <- data_process_norfolk(sl.data.path=file.path(data.path, 'tide_gauge_SewellsPoint'),
                                   temp.data.path=file.path(data.path, 'temperature'))
  saveRDS(dat_list, file.path(output.path, filename.data))
} else {
  dat_list <- readRDS(file.path(output.path, filename.data))
}

# prune temperature data to match tide gauge data
dat_calib <- merge(dat_list$gev, dat_list$temps$forcing)

# set up parameter name list
if (model == 'stationary') {
  parnames <- c('lambda', 'sigma', 'xi')
} else if (model == 'nonstationary') {
  parnames <- c('lambda0', 'lambda1', 'sigma', 'xi')
}

# if method is MLE, compute maximum likelihood estimates and save. This must be done
# before MCMC.
if (method == 'MLE') {
  filename.mle <- paste('norfolk_mle-',model,'.rds', sep='')
  
  if (file.exists(file.path(output.path, filename.mle))) {
    print('MLE Parameters already found!')
  } else {
    print('starting DE optimization for MLE parameters for Norfolk...')
  
    if (model == 'stationary') {
      aux_data <- NULL
    } else if (model == 'nonstationary') {
      aux_data <- dat_calib$temp
    }
    deoptim_out <- DEoptim(fn=neg_log_likelihood, 
                           lower=bounds.lower,
                           upper=bounds.upper,
                           DEoptim.control(NP=NP.deoptim, itermax=niter.deoptim, F=F.deoptim, CR=CR.deoptim, trace=FALSE),
                           parnames=parnames, model=model, dat=dat_calib$ann.max, aux_data=aux_data)
  
    est.ml <- deoptim_out$optim$bestmem
    names(est.ml) <- parnames
    
    saveRDS(est.ml, file.path(output.path, filename.mle))
  }
} else {
  # otherwise, load the priors if available, or fit them and save if not
  filename.prior <- paste('tide_gauge_priors-', model, '.rds', sep='')
  
  if (file.exists(file.path(output.path, filename.prior))) {
    print('Priors already found! Loading...')
    priors <- readRDS(file.path(output.path, filename.prior))
  } else {
    print('fitting prior distributions...')
    source(file.path(R.path, 'fit_priors.R'))
    priors <- fit_priors_stations(tide.data.path=file.path(data.path, 'tide_gauge_long'),
                                  temp.data.path=file.path(data.path, 'temperature'),
                                  model=model)
    saveRDS(priors, filename.prior)
  }
  
  # a preliminary run (mode == 'prelim') is required to get an estimate of the jump covariance matrix
  # as well as initial points for the production MCMC chains
  if (method == 'prelim') {
    # load MLE fit
    filename.mle <- paste('norfolk_mle-',model,'.rds', sep='')
    if (!file.exists(file.path(output.path, filename.mle))) {
      stop("need to run in MLE mode first!")
    } else {
      mle.fit <- readRDS(file.path(output.path, filename.mle))
    }
    
    rate.accept_many <- 0.234 # optimal acceptance rate for multiple parameters
    rate.accept_one <- 0.44 # optimal acceptance rate for a single parameter
    rate.accept <- rate.accept_many + (rate.accept_one - rate.accept_many)/length(parnames)
    
    # set step size to be 5% of the standard deviation of the prior distribution
    stepsize <- numeric(length(parnames))
    for (par in parnames) {
      if (priors[[par]]$type == 'normal') {
        stepsize[match(par, parnames)] <- 0.05*priors[[par]]$sd
      } else if (priors[[par]]$type == 'gamma') {
        stepsize[match(par, parnames)] <- 0.05*sqrt(priors[[par]]$shape*priors[[par]]$scale^2)
      }
    }
    adapt.start <- max(500, round(0.05*n.iter))
    
    print(paste('starting preliminary calibration for', model, 'model...', sep=' '))
    
    # start MCMC from the MLE values
    if (model == 'stationary') {
      aux_data <- NULL
    } else if (model == 'nonstationary') {
      aux_data <- dat_calib$temp
    }

    amcmc_prelim <- MCMC(log_posterior_gev, n.iter, init=mle.fit, scale=stepsize, gamma=0.67, 
                         list=TRUE, n.start=adapt.start, acc.rate=rate.accept, parnames=parnames, priors=priors,
                         model=model, dat=dat_calib$ann.max, aux_data=aux_data)
    
    # save preliminary MCMC results
    saveRDS(amcmc_prelim, file.path(output.path, paste('norfolk_MCMC_prelim-', model, '.rds', sep='')))
   print('...done!') 
  } else {
    # read in preliminary MCMC results
    amcmc_prelim <- readRDS(file.path(output.path, paste('norfolk_MCMC_prelim-', model, '.rds', sep='')))
    
    init.value <- amcmc_prelim$samples[nrow(amcmc_prelim$samples),]
    
    adapt.start <- max(500, round(0.05*n.iter))
    if (model == 'stationary') {
      aux_data <- NULL
    } else if (model == 'nonstationary') {
      aux_data <- dat_calib$temp
    }
    
    rate.accept_many <- 0.234 # optimal acceptance rate for multiple parameters
    rate.accept_one <- 0.44 # optimal acceptance rate for a single parameter
    rate.accept <- rate.accept_many + (rate.accept_one - rate.accept_many)/length(parnames)
    
    stepsize <- amcmc_prelim$cov.jump
    
    amcmc_prod <- MCMC.parallel(log_posterior_gev, n.iter, init=init.value, n.chain = n.chain, n.cpu = n.cpu, 
                                scale=stepsize, gamma=0.67, list=TRUE, n.start=adapt.start, acc.rate=rate.accept, 
                                packages=c('extRemes'),
                                parnames=parnames, priors=priors, model=model, dat=dat_calib$ann.max, aux_data=aux_data)
    
    saveRDS(amcmc_prod, file.path(output.path, paste('norfolk_MCMC-', model, '.rds', sep='')))
  }
}