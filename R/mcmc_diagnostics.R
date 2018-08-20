###################################################################################################
# mcmc_diagnostics.R: this file contains functions for MCMC convergence diagnostics and		  #
# 		      requires the coda package.						  #
###################################################################################################
library(coda)

# Gelman-Rubin diagnostics
# assume that you pass a list of MCMC output objects, with each row corresponding to an iteration.
# returns the multivariate gr diagnostic for every 1/20th of each chain length
gr_diag <- function(chains) {
  if (is.list(chains)) {
    iter_length <- nrow(chains[[1]]$samples)
    test_seq <- seq(1, iter_length, by = 0.05*iter_length)
    sapply(test_seq, function(m) as.numeric(gelman.diag(lapply(chains, function(l) as.mcmc(l$samples[1:m,])))[[2]]))
  } else {
    iter_length <- nrow(chains$samples)
    test_seq <- seq(1, iter_length, by = 0.05*iter_length)
    sapply(test_seq, function(m) as.numeric(gelman.diag(as.mcmc(chains$samples[1:m,]))[[2]]))
  }
}

# monitor 5%, 50%, 95% quantiles of each variabile for drifting
quant_diag <- function(mcmc_list) {
  iter_length <- nrow(chains[[1]])
  test_seq <- seq(1, iter_length, by = 0.05*iter_length)
  lapply(mcmc_list, function(l) apply(l$samples, 2, quantile, probs=c(0.05, 0.5, 0.95)))
}
  