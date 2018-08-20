# SPSLAM (Sewell's Point Sea Level Assessment Model)
Sea-level rise and storm surge analysis and modeling for Sewell's Point, VA.

## Overview
This repository contains code to fit models for sea-level rise and storm surge at Sewell's Point, VA, as well as to replicate the data analysis and figures from Srikrishnan et al (2018), "Investing in Science and Using the Results to Improve Climate Risk Management." The code assumes that MCMC chains will be generated in parallel, and some tweaking is required to run the code in serial.

## Requirements
### Data
All required data for the linear sea-level rise scenario and the storm surge models are available in the data/ directory. For the non-linear sea-level rise scenarios, BRICK (https://github.com/scrim-network/BRICK) and its calibration file are required. 

### R packages
The following R packages are required:
* DEoptim (finding maximum likelihood estimates)
* ncdf4 (data I/O)
* coda (for the Gelman-Rubin diagnostic)
* adaptMCMC (running adaptive MCMC)
* parallel (generating multiple MCMC chains in parallel)
* extRemes (likelihood and sampling functions for GEV models)
* plyr (data frame manipulation)
* reshape2 (melting data frames)
* ggplot2 (plotting)
* scales (additional scales for axes)
* gridExtra (additional layout functions for ggplot2)
* ggpubr (additional ggplot2 themes)

## Workflow
After installing BRICK, downloading its calibration file, and installing the required R packages, reproduction involves the following steps:
1. Process data and calibrate storm surge models: 
`Rscript R/calibrate_model.R`
2. Calibrate linear sea-level rise model:
`Rscript R/linear_slr_MCMC.R`
3. Generate ensemble of sea-level rise futures and compute 100-year return levels:
`Rscript R/compute_rl.R`
4. Plot figures:
`Rscript R/plot_figure1.R`
`Rscript R/plot_figure2.R`
