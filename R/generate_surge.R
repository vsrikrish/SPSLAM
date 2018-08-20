###################################################################
# generate_surge.R                                                #
#   This function generates an ensemble of annual storm surge     #
#     values at Sewell's Point, VA. The surges are assumed        #
#     to follow a GEV distribution with a non-stationary          #
#     location parameter, which depends on the global mean        #
#     temperature anomaly relative to 1901--2000.                 #
#                                                                 #
#   The main function to be called is rand_surge(n, temps),       #
#     which returns a list or matrix of n sampled surge values    #
#     for each year in the temperature anomaly dataset.           #
#                                                                 #
#     It assumes that the temperatures are passed either as a     #
#     data frame or a list of data frames (if multiple            #
#     temperature scenarios are used). Each data frame should     #
#     have two columns: the first column should contain the       #
#     years, and the second the temperature anomaly. The years    #
#     are used as the resulting matrix column names. If a list    #
#     of multiple temperature scenarios is provided, the output   #
#     is a list of matrices.                                      #
#                                                                 #
#   This function requires the 'SewellsPoint-MCMC.rds'            #
#     file, which is assumed to be placed in an 'output/' subdir. #
#     It also requires the extRemes package.                      #
###################################################################


# load extRemes package
library(extRemes)

###################################################################
# rand_surge_singleSOW(n, temps, params):                         #
#   This function samples surge values for each individual        #
#     temperature scenario.                                       #
#   Inputs:                                                       #
#     n (numeric): number of samples for each year                #
#     temps (data frame): 2-column data frame where the           #
#         years are the first column and the temperature          #
#         anomalies are the second column.                        #
#     params: matrix of MCMC parameter samples.                   #
#   Output:                                                       #
#     n x num_years matrix, where num_years is the number of      #
#       rows of the temperature series.                           #
###################################################################
rand_surge_singleSOW <- function(n, temps, params, stationary=T) {
  # generate parameter ensembles
  
  if (!stationary) {
    lambda <- vapply(temps[,2],
                FUN=function(t) {
                  params[, match('lambda0', colnames(params))] +
                  params[, match('lambda1', colnames(params))] * t},
                FUN.VALUE=numeric(nrow(params)))
  } else {
    lambda <- params[, match('lambda', colnames(params))]
  }
    sigma <- params[, match('sigma', colnames(params))]
    xi <- params[, match('xi', colnames(params))]
    
    # generate a matrix of surge samples; columns are by year
    surges <- vapply(1:length(temps[,2]),
                FUN=function(i) {
                  revd(n, loc=lambda[,i], scale=sigma, shape=xi, type='GEV')},
                FUN.VALUE=numeric(nrow(params)))

  # name the columns of the surges using years
  colnames(surges) <- temps[,1]
  
  # return surge matrix
  surges
}

###################################################################
# rand_surge(n, temps):                                           #
#   This function samples surge values for each element of a      #
#     list of temperature scenarios. If only temperature          #
#     scenario is provided, it returns a single matrix instead    #
#     of a list.                                                  #
#   Inputs:                                                       #
#     n (numeric): number of samples for each year.               #
#     temps (data frame or list of data frames):                  #
#         2-column data frame or list of 2-column data frames,    #
#         where each data frame has the years as the first        #
#         column and the temperature anomalies as the second      #
#         column.                                                 #
#   Output:                                                       #
#     list of n x num_years matrices, where num_years is the      #
#       number of rows of each temperature series. If a single    #
#       temperature data frame is provided, a single matrix       #
#       is returned.                                              #
###################################################################
rand_surge <- function(n, temps, stationary=T) {
  ## load MCMC output
  # this assumes the MCMC file "SewellsPoint-MCMC.rds" is located in a
  # "data" subdirectory and that the function is called from the main
  # working directory
  if (!stationary) {
    params <- readRDS('data/SewellsPoint-MCMC.rds')[['nonstationary']]
  } else {
    params <- readRDS('data/SewellsPoint-MCMC.rds')[['stationary']]
  }

  ## sample ensemble members
  idx <- sample(1:nrow(params), n, replace = T)
  
  ## sample over each temperature list if multiple temperature scenarios
  #   are provided
  if (is.data.frame(temps)) {
    rand_surge_singleSOW(n, temps, params[idx,], stationary)
  } else {
    lapply(temps,
            FUN=function(l) {
              rand_surge_singleSOW(n, l, params[idx,], stationary)
            })
  }
}
