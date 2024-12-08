## ---------------------------
##
## Script name: 06_multispeciesoccupancy.R 
##
## Purpose of script: Basic multi-species occupancy models in R with JAGS
##
## Author: James E Paterson
##
## Date Created: 2024-12-08
##
## ---------------------------


## ----loaddata--------------------------------------------------
# Load libraries
library(dplyr) # for data organization
library(ggplot2) # for plotting

# Load arrary of community detection history 
# 100 x (sites) x 10 columns (surveys) x 18 species
# It is the same structure as a single species detection history, 
# but stacked with 17 additional species
load(file = "community_detection_history.RData")

# Structure
str(community_detection_history)

# Look at part of one species' detection history (species 10)
head(community_detection_history[,,10])

# Restructure detection history data to summarize how many surveys a species was detected at
y <- matrix(data = NA, 
            nrow = nrow(community_detection_history), 
            ncol = dim(community_detection_history)[3])

# Use a simple loop to fill in y
# for each site row
for(i in 1:nrow(community_detection_history)){ 
  # for each species column
  for(k in 1:dim(community_detection_history)[3]){ 
    # sum the number of times a species was detected
    y[i,k] <- sum(community_detection_history[i,,k]) 
  }
}



## ----simplemsom------------------------------------------------
# Load R to jags library
# If installing for the first time, you also need to install jags.
# https://mcmc-jags.sourceforge.io/
library(jagsUI)

# Convert detection history to numeric
z <- (y > 0)*1  # *1 converts to numeric, keeps matrix
z[z == 0] <- NA

# Pack up data for JAGS
msom_jags_data <- list(y = y, # y = the detection history matrix
                 nSobs = ncol(y), # nSobs = number of species observed
                 nSites = nrow(y),  # nSites = number of sites
                 nOcc = ncol(community_detection_history), # nOcc = number of surveys
                 z = z) # z =  input detection data

# Look at structure
str(msom_jags_data)

# Hierarchical model
# For fitting basic multi-species occupancy models, we can think of the components as
# 1. Ecological model
# 2. Observation model
# 3. Priors that describe species level variation in psi and p
# 4. Hyperpriors that describe parameters for the community 

# Create and save file name "msom_simple.jags"
# Identify filepath of model file
msom_simple <- tempfile()

#Write model to file
writeLines("
model{
  for(k in 1:nSobs){  # Loop through species
    # Likelihood
    for(i in 1:nSites) { # Loop through sites
      # Ecological model
      z[i, k] ~ dbern(psi[k])
      # Observation model
      y[i, k] ~ dbin(p[k] * z[i, k], nOcc)
    }
    
    # Priors for Psi (species level)
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    psi[k] <- ilogit(lpsi[k])
    
    # Priors for p (species level)
    lp[k] ~ dnorm(mu.lp, tau.lp)
    p[k] <- ilogit(lp[k])
  }
  
  # Hyperpriors for Psi (community level)
  psi.mean ~ dbeta(1, 1)
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- 1/sd.lpsi^2
  
  # Hyperpriors for p (community level)
  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- 1/sd.lp^2
}
", con = msom_simple)

# Vector of monitored parameters
wanted <- c("psi.mean", "p.mean", "psi", "p")
# Could add other parameters in model: "mu.lpsi", "sd.lpsi", "tau.lpsi",
# "mu.lp", "sd.lp", "tau.lp"

# model run takes <1 min on a macbook
msom_simple_out <- jags(msom_jags_data, NULL, wanted, msom_simple,
                        n.chains = 3, n.adapt = 100, n.iter = 10000, 
                        n.burnin = 5000, n.thin = 2, parallel = TRUE, DIC = FALSE)


## ----checkdiagnostics ------------------------------------------
# Plot traces and density estimates 
# It will show many pages of plots
 plot(msom_simple_out)

## ----parameterposteriordist----------------------------------------------------------------------
# Summary of results (using head() to just show higher levels)
# Full summary includes psi and p estimates for each species
head(summary(msom_simple_out))

## ----plotpsiandp-------------------------------------------------------------------------------------
# Put psi estimates in a dataframe
psi_species_estimates <- summary(msom_simple_out) %>%
  as.data.frame(.) %>%
  mutate(parameter = row.names(.)) %>%
  # Filtering to only the estimates of psi for each species 
  # species psi parameters in our model take the form "psi[1]"
  filter(parameter %in% c(paste0("psi","[",c(1:18),"]"))) %>%
  arrange(mean) %>%
  mutate(plot_order = 1:nrow(.))
  
# Plot psi estimates
ggplot(psi_species_estimates, aes(x = as.factor(plot_order), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
  labs(x = "Species", y = "Occupancy estimate") +
  scale_x_discrete(labels = psi_species_estimates$parameter, 
                   breaks = order)+ 
  scale_y_continuous(limits = c(0, 1.0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  coord_flip() +
  theme_classic()
  
# Put p estimates in a dataframe
p_species_estimates <- summary(msom_simple_out) %>%
  as.data.frame(.) %>%
  mutate(parameter = row.names(.)) %>%
  # Filtering to only the estimates of psi for each species 
  # species psi parameters in our model take the form "psi[1]"
  filter(parameter %in% c(paste0("p","[",c(1:18),"]"))) %>%
  arrange(mean) %>%
  mutate(plot_order = 1:nrow(.))
  
# Plot p estimates
ggplot(p_species_estimates, aes(x = as.factor(plot_order), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
  labs(x = "Species", y = "Detection probability estimate") +
  scale_x_discrete(labels = p_species_estimates$parameter, 
                   breaks = order)+ 
  scale_y_continuous(limits = c(0, 1.0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  coord_flip() +
  theme_classic()
