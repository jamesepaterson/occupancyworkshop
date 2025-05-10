## ---------------------------
##
## Script name: 08_msomundetectedspecies.R 
##
## Purpose of script: Multi-species occupancy models in R: covariates and estimating species richness
##
## Author: James E Paterson
##
## Date Created: 2025-05-09
##

## ----loaddata-----------------------------------------------------

# Load libraries
library(dplyr) # for data organization
library(ggplot2) # for plotting

# Load arrary of community detection history 
# 100 x rows (sites) x 18 columns (species) with values for the number of detections in 10 surveys
y <- read.csv("community_detection_3visits.csv") %>%
  dplyr::select(-X)

# Convert detection history to numeric
z <- (y > 0)*1  # *1 converts to numeric, keeps matrix
z[z == 0] <- NA

# Load sitecov, effort, observers, 
site.cov <- read.csv(file = "site_cov.csv") %>%
  dplyr::select(-X)


## ----msomCAugmented-----------------------------------------------

library(jagsUI)

naug <- 50 - ncol(y) # 18 observed species, we'll augment the data with an additional 32 species for a total of 50

# Updated y; matrix of detection histories that includes detected and non-detected species
yaug <- cbind(y, matrix(0, nrow(y), naug))

# Initial w vector (1 for detected species, NA for augmented undetected species)
w <- c(rep(1, ncol(y)), rep(NA, naug))

# Initial values for latent state occupancy (1 for detected, otherwise NA)
z <- (yaug > 0)*1
z[z == 0] <- NA

# Package data for jags
msom_augmented_jags_data <- list(y = yaug,# y = the detection history matrix
                                  nSobs = ncol(yaug), # The total number of species observed + augmented
                                  nSites = nrow(yaug), # The total number of sites
                                  nOcc = 3, # Each site visited 10 times
                                  z = z, # z = latent state occupancy
                                  w = w, # w = binary identifier of whether a species is part of the community
                                   # Covariates (already standardized)
                                  forest = site.cov$forest,
                                  agri = site.cov$agri)

# Hierarchical model
# For fitting basic multi-species occupancy models, we can think of the components as
# 1. an ecological model
# 2. an observation model
# 3. Priors that describe species level variation in psi and p
# 4. hyperpriors that describe parameters for the community 
# Create a temporary file "msom_augmented"
# Identify filepath of model file
msom_augmented <- tempfile()

# Write model to file
writeLines("
model{
  for(k in 1:nSobs){  # Loop through species, including augmented data
  w[k] ~ dbern(omega)
    # Likelihood
    for(i in 1:nSites) {
       # Ecological model
      logit(psi[i, k]) <- b0[k] + bforest[k]*forest[i] + bagri[k]*agri[i]
       z[i, k] ~ dbern(psi[i, k]*w[k]) # Addition of w[k] qualifying whether species are at a site
       # Observation model
       y[i, k] ~ dbin(p[k] * z[i, k], nOcc)
    }

    # Priors (species level)
    b0[k] ~ dnorm(mu.b0, tau.b0)
    # mu.eta[k] <- mu.lp + rho * sd.lp/sd.b0 *
    #     (b0[k] - mu.b0)
    # lp[k] ~ dnorm(mu.eta[k], tau.eta)
    lp[k] ~ dnorm(mu.lp, tau.lp)
    p[k] <- ilogit(lp[k])

    bforest[k] ~ dnorm(mu.forest, tau.forest)
    bagri[k] ~ dnorm(mu.agri, tau.agri)

  }

  # Hyperpriors (community level)
  omega ~ dbeta(0.001, 1) # Link (2013) suggested strong prior Beta(0.001, 1)
  
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
  
  # Derived species richness
  N <- sum(w)
}
", con = msom_augmented)

# Parameters monitored (greatly reducing to focus on omega, N, z, w)
wanted <- c("omega", "N", 
            "z", "w")

# Run model
msom_augmented_out <- jags(msom_augmented_jags_data, 
                           NULL, wanted, msom_augmented,
                           n.chains = 3, n.iter = 25000, 
                           n.burnin = 5000, n.thin = 5, 
                           parallel = TRUE, DIC = FALSE)

# Summary
modelSummary <- summary(msom_augmented_out)

## ----omegaNresults----------------------------------------------------------------------------------

# Posterior distribution of Omega
modelSummary["omega",]

# Posterior distribution of N
modelSummary["N",]

# Plot distribution of N
ggplot(data = data.frame(N = msom_augmented_out$sims.list$N), aes(x = N)) +
  # geom_histogram(fill = "white", col = "black") +
  geom_density(col = "dark green", adjust = 4, fill = "light green") +
  geom_vline(aes(xintercept = median(msom_augmented_out$sims.list$N)),
            color = "dark green", linetype = "dashed", size = 1)+
  geom_vline(aes(xintercept = 18),
            color = "black", linetype = "dashed", size = 1)+
  scale_x_continuous(breaks = c(18:30), limits = c(18, 30))+
  labs(x = "Estimated number of species in community (detected and undetected)",
       y = "Density") +
  theme_classic()


## ----speciesRichnessEstimates-----------------------------------------------------------------------

# Species occupancy state at each site
z_results <- msom_augmented_out$sims.list$z

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
# Horizontal line at total number of unique species observed (18)
ggplot(richnessDF, aes(x = forest, y = meanRichness)) +
  geom_point(pch = 16) +
  geom_errorbar(aes(ymin = lowerRichness, ymax = upperRichness), width = 0) +
  labs(y = "Predicted species richness", x = "Forest amount (scaled)") +
  geom_hline(aes(yintercept = 18), linetype = "dashed") +
  theme_classic()

## END