# Workshop: occupancy models in R

A workshop on occupancy models in R. This workshop uses the 'unmarked', 'RPresence', and 'jagsUI' packages to construct occupancy models in R. All examples use included simulated datasets.

This workshop covers:

* single season occupancy models and goodness-of-fit tests
* building and comparing multiple models
* model-averaging predictions for occupancy models
* fitting dynamic (multi-season) occupancy models
* fitting implicit dynamics occupancy models
* fitting multi-species occupancy models with JAGS

Contents:

* 01_installingoccupancypackages.R checks if required packages are installled. Installs them if they are not.
* 02_singleseasonoccupancy.R contains code for introduction to occupancy models
* 03_singleseasonoccupancy_part2_modelcomparison.R contains code for building and comparing single season occupancy models
* 04_dynamicoccupancy.R contains code for building dynamic (multi-season) occupancy models
* 05_implicitdynamicsoccupancy.R contains code for building implicit dynamic occupancy models
* 06_multispeciesoccupancy.R contains code for building basic multi-species occupancy models with JAGS
* 07_multispeciesoccupancycovariates.R contains code for multi-species occupancy models with covariates and estimating species richness
* detection_history.csv contains simulated detection history data
* community_detection_history.RData contains simulated community detection data for 18 species
* community_detection_history_covariates.RData contains simulated community detection data for 18 species for covariate tutorial
* site_cov.csv contains simulated site covariates (forest cover and agriculture cover)
* effort.csv contains simulated observation covariate for search effort per survey
* observers.csv contains simulated observation covariate for number of observers per survey
* dynamic_detection_history.csv contains simulated detection history data used in 04_dynamicoccupancy
* dynamic_site_cov.csv contains simulated site covariates (wetland, temp, prec) used in 04_dynamicoccupancy
* dynamic_effort.csv contains simulated observation covariate for search effort per survey used in 04_dynamicoccupancy
