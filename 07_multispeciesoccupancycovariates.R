## ---------------------------
##
## Script name: 07_multispeciesoccupancycovariates.R 
##
## Purpose of script: Multi-species occupancy models in R: covariates and estimating species richness
##
## Author: James E Paterson
##
## Date Created: 2025-01-19
##

## ----loaddata --------------------------------
# Load libraries
library(dplyr) # for data organization
library(ggplot2) # for plotting
library(tidyr) # for pivoting data

# Load arrary of community detection history 
# 100 x (sites) x 10 columns (surveys) x 18 species
# It is the same structure as a single species detection history, 
# but stacked with 17 additional species
load(file = "community_detection_history_covariates.RData")

community_detection_history <- community_detection_history_covariates
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

# Convert detection history to numeric
z <- (y > 0)*1  # *1 converts to numeric, keeps matrix
z[z == 0] <- NA

# Load sitecov, effort, observers, 
site.cov <- read.csv(file = "site_cov.csv") %>%
  dplyr::select(-X)

## ----msomCovariates--------------------------
library(jagsUI)

# Package data for jags, which must specify covariates
msom_covariates_jags_data <- list(y = y,# y = the detection history matrix
                                  nSobs = ncol(y), # The total number of species observed
                                  nSites = nrow(y),  # The total number of sites along the route
                                  nOcc = ncol(community_detection_history), # Each site visited 10 times
                                  z = z,
                                   # Covariates (already standardized)
                                  forest = site.cov$forest,
                                  agri = site.cov$agri)

# Hierarchical model
# For fitting basic multi-species occupancy models, we can think of the components as
# 1. an ecological model
# 2. an observation model
# 3. Priors that describe species level variation in psi and p
# 4. hyperpriors that describe parameters for the community 
# Create a temporary file "msom_covariates"
#Identify filepath of model file
msom_covariates <- tempfile()

#Write model to file
writeLines("
model{
  for(k in 1:nSobs){  # Loop through species
    # Likelihood
    for(i in 1:nSites) {
       # Ecological model
      logit(psi[i, k]) <- b0[k] + bforest[k]*forest[i] + bagri[k]*agri[i]
       z[i, k] ~ dbern(psi[i, k])
       # Observation model
       y[i, k] ~ dbin(p[k] * z[i, k], nOcc)
    }

    # Priors (species level)
    b0[k] ~ dnorm(mu.b0, tau.b0)
    mu.eta[k] <- mu.lp + rho * sd.lp/sd.b0 *
        (b0[k] - mu.b0)
    lp[k] ~ dnorm(mu.eta[k], tau.eta)
    p[k] <- ilogit(lp[k])

    bforest[k] ~ dnorm(mu.forest, tau.forest)
    bagri[k] ~ dnorm(mu.agri, tau.agri)

  }

  # Hyperpriors (community level)
  
  b0.mean ~ dbeta(1, 1)
  mu.b0 <- logit(b0.mean)
  sd.b0 ~ dunif(0, 5)
  tau.b0 <- 1/sd.b0^2

  mu.forest ~ dunif(-5, 5)
  sd.forest ~ dunif(0, 5)
  tau.forest <- 1/sd.forest^2
  
  mu.agri ~ dunif(-5, 5)
  sd.agri ~ dunif(0, 5)
  tau.agri <- 1/sd.agri^2

  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- 1/sd.lp^2

  rho ~ dunif(-1, 1)
  tau.eta <- tau.lp/(1 - rho^2)
  
}
", con = msom_covariates)

# Parameters monitored
wanted <- c("b0.mean", "mu.b0", "sd.b0", 
            "psi", "b0", "bforest",
            "mu.forest", "sd.forest",
            "mu.agri", "sd.agri",
            "p.mean", "sd.lp", "z")

# model run takes under 2 min
system.time(
msom_covariates_out <- jags(msom_covariates_jags_data, NULL, wanted, msom_covariates,
             n.chains = 3, n.iter = 10000, 
             n.burnin = 1000, n.thin = 5, parallel = TRUE, DIC = FALSE)
)

# Summary
modelSummary <- summary(msom_covariates_out)

## ----checkdiagnostics------------------------
# Plot traces and density estimates 
# (not showing because it will show many pages of plots)
# plot(msom_covariates_out)

# Summary
head(summary(msom_covariates_out))

## ----plotcovariateeffects--------------------
forestCoefficients <- summary(msom_covariates_out) %>%
  as.data.frame(.) %>%
  mutate(parameter = row.names(.)) %>%
  # Filtering to only the estimates of psi for each species 
  # species psi parameters in our model take the form "psi[1]"
  filter(parameter %in% c(paste0("bforest","[",c(1:18),"]"))) %>%
  arrange(mean) %>%
  mutate(plot_order = 1:nrow(.))

# Plot forest coefficient estimates
ggplot(forestCoefficients, aes(x = as.factor(plot_order), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = `25%`, ymax = `75%`), width = 0, size = 2) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
  labs(x = "Species", y = "Forest Coeffiecient Estimate on Occupancy") +
  scale_x_discrete(labels = forestCoefficients$parameter, 
                   breaks = order)+ 
  scale_y_continuous(limits = c(-0.25, 1.15), breaks = c(-0.25, 0, 0.25, 0.5, 1.0)) +
  coord_flip() +
  theme_classic()

# Extract mean predicted relationships

# Create a sequence to predict on from the minimum to maximum of a covariate
forestSequence <- seq(from = range(site.cov$forest)[1], 
                      to = range(site.cov$forest)[2], 
                      length.out = 101)
  
# Get the predictions for each value of forestSequence
predForest <- matrix(NA, 101, ncol(y))
for(i in 1:101){
  predForest[i,] <- plogis(msom_covariates_out$mean$b0 +
      msom_covariates_out$mean$bforest * forestSequence[i])
}

# Add forestSequence then pivot longer
predForest<- predForest %>%
  as.data.frame(.) %>%
  mutate(forest = forestSequence) %>%
  pivot_longer(., cols = 'V1':'V18', names_to = "species", values_to = "occupancy")

# Plot relationship between species occupaancy and forest
ggplot(predForest, aes(x = forest, y = occupancy, group = species)) +
    geom_path(col = "black") +
    labs(x = "Forest amount (scaled)", y = "Predicted occupancy") +
    theme_classic()

## ----covariatesRichness--------------------------------------------------------
# Species occupancy state at each site
z_results <- msom_covariates_out$sims.list$z

# Calculate average species richness per site
ndraws <- dim(z_results)[1] # Sample size (number of draws from posterior distribution)
richnessMatrix <- matrix(nrow = 100, ncol = ndraws) # Matrix to hold results

# For each draw from the posterior distribution
for(i in 1:ndraws){
  # Calculate the species richness by summing species
  # Each row is a site, each column is one draw
  richnessMatrix[,i] <- rowSums(z_results[i,,])
}

richnessDF <- richnessMatrix %>%
  as.data.frame(.) %>%
  # Calculate the average, SD, lower (2.5%) and upper (95.5%) of the distribution per site
  mutate(meanRichness = rowMeans(dplyr::select(., starts_with("V")), na.rm = TRUE),
         sdRichness = apply(X = richnessMatrix, MARGIN = 1, FUN = sd),
         lowerRichness = apply(X = richnessMatrix, MARGIN = 1, 
                               FUN = quantile, probs = 0.025),
         upperRichness = apply(X = richnessMatrix, MARGIN = 1, 
                               FUN = quantile, probs = 0.975),
         sites = as.character(1:100)) %>%
  # Join to the site dataframe to have the covariate
  left_join(., site.cov %>%
              mutate(sites = as.character(sites)),
            by = "sites") %>%
  # Limit the data frame to just the variables we need
  dplyr::select(sites, meanRichness, sdRichness,
                lowerRichness, upperRichness,
                forest, agri)

# Plot richness vs forest
ggplot(richnessDF, aes(x = forest, y = meanRichness)) +
  geom_point(pch = 16) +
  geom_errorbar(aes(ymin = lowerRichness, ymax = upperRichness), width = 0) +
  labs(y = "Predicted species richness", x = "Forest amount (scaled)") +
  theme_classic()