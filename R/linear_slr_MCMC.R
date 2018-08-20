######################################################################################
# linear_slr_MCMC.R: file to do MCMC for the linear trend at Sewell's Point	     #
# 		     This file fits a linear regression model for the mean SLR	     #
#		     to use only the historical data, ignoring structural  	     #
#		     contributors such as ice sheet dynamics.			     #
		     It requires the adaptMCMC package.				     #
######################################################################################

library(adaptMCMC)

# parameters
n.iter <- 5e5

# likelihood for basic linear regression with normal coefficient priors and gamma variance priors
log_lik_regress <- function(params, parnames, dat, cov) {
  beta0 <- params[match('beta0', parnames)]
  beta1 <- params[match('beta1', parnames)]
  sigma.sq <- params[match('sigma', parnames)]^2

  -0.5*length(dat)*log(2*pi*sigma.sq) - 0.5/sigma.sq*(t(dat-(beta0+beta1*cov))%*%(dat-(beta0+beta1*cov)))
}

# priors for the coefficients are centered at the MLE and are relatively diffuse (depending on scale)
# gamma(2,2) prior for the sd is also based on the MLE fit
log_prior_regress <- function(params, parnames) {
  beta0 <- params[match('beta0', parnames)]
  beta1 <- params[match('beta1', parnames)]
  sigma <- params[match('sigma', parnames)]
  
  dnorm(beta0, mean=-288, sd=100, log=TRUE) + dnorm(beta1, mean=4.6, sd=1, log=TRUE) + dgamma(sigma, shape=2, scale=10, log=TRUE)
}

log_post_regress <- function(params, parnames, dat, cov) {
  log.pri <- log_prior_regress(params, parnames)
  if (is.finite(log.pri)) {
    log.lik <- log_lik_regress(params, parnames, dat, cov)
  } else {
    log.lik <- 0
  }
  log.pri + log.lik
}

dat <- readRDS('output/processed_norfolk.rds')
dat[['lin']]$yridx <- dat[['lin']]$year - dat[['lin']]$year[1]
amcmc_out <- MCMC.parallel(log_post_regress, n.iter, n.chain=10, n.cpu=parallel::detectCores(),
                           init=c(-288, 4.6, 30), acc.rate=0.37, list=TRUE, n.start=0.05*n.iter,
                           parnames=c('beta0', 'beta1', 'sigma'), dat=dat[['lin']]$mean, cov=dat[['lin']]$yridx)
saveRDS(amcmc_out, 'output/MCMC-norfolk_linear_slr.rds')