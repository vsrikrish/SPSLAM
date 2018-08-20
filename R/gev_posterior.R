# functions to compute log likelihoods for GEV fit

library(extRemes)

# compute log prior
# priors should be passed as a list with each entry as follows:
#   parameter: name of parameter;
#   other parameters named appropriately, such as "mean" and "sd" for normal family priors
log_prior <- function(params,   # vector of sampled parameters
                      parnames, # character vector of parameter names
                      priors,   # list of prior settings
                      model     # 'stationary' or 'nonstationary
){
  prior.log <- 0  # initialize prior log-likelihood to 0
  # for each parameter, compute likelihood
  for (par in parnames) {
    par.val <- params[match(par, parnames)]
    if (priors[[par]]$type == 'normal') {
      prior.log <- prior.log + dnorm(par.val, mean=priors[[par]]$mean, sd=priors[[par]]$sd, log=TRUE)
    } else if (priors[[par]]$type == 'lognormal') {
      prior.log <- prior.log + dlnorm(par.val, meanlog=priors[[par]]$meanlog, sdlog=priors[[par]]$sdlog, log=TRUE)
    } else if (priors[[par]]$type == 'gamma') {
      prior.log <- prior.log + dgamma(par.val, shape=priors[[par]]$shape, scale=priors[[par]]$scale, log=TRUE)
    }
  }
  # return
  prior.log
}

log_likelihood_gev <- function(params,   # vector of sampled parameters
                               parnames, # character vector of parameter names
                               model,    # 'stationary' or 'nonstationary'
                               dat,      # vector of block maxima
                               aux_data = NULL # covariates for non-stationary processes 
){
  if (model == 'stationary') {
    lambda <- params[match('lambda', parnames)]
    sigma <- params[match('sigma', parnames)]
    xi <- params[match('xi', parnames)]
  } else if (model == 'nonstationary') {
    lambda <- params[match('lambda0', parnames)] + params[match('lambda1', parnames)]*aux_data
    sigma <- params[match('sigma', parnames)]
    xi <- params[match('xi', parnames)]
  }
  # compute log likelihood
  sum(devd(dat, loc=lambda, scale=sigma, shape=xi, log=TRUE, type='GEV'))
}
  
log_posterior_gev <- function(params,   # vector of sampled parameters
                              parnames, # character vector of parameter names
                              priors,   # list of prior settings
                              model,    # 'stationary' or 'nonstationary'
                              dat,      # vector of block maxima
                              aux_data = NULL # covariates for non-stationary processes 
){
  # throw error message if covariates are missing for nonstationary model
  if ((model == 'nonstationary') & is.null(aux_data)) {
    stop("aux_data needs to exist for nonstationary models!")
  }
  # compute prior density for parameters
  prior.log <- log_prior(params, parnames, priors, model)
  
  # if parameters are plausible wrt priors, compute likelihood
  if (is.finite(prior.log)) {
    # get parameter values
    lik.log <- log_likelihood_gev(params, parnames, model, dat, aux_data)
  } else {
    # if parameters are not plausible, just set log likelihood to zero
    lik.log <- 0
  }
  
  # return the log posterior is the sum of the log prior and the log likelihood
  prior.log + lik.log
  
}

# return the negative log likelihood for computing the MLE
neg_log_likelihood <- function(params,   # vector of sampled parameters
                               parnames, # character vector of parameter names
                               model,    # 'stationary' or 'nonstationary'
                               dat,      # vector of block maxima
                               aux_data = NULL # covariates for non-stationary processes 
){
  -1*log_likelihood_gev(params, parnames, model, dat, aux_data)
}