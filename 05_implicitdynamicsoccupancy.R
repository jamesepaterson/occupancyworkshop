## ---------------------------
##
## Script name: 05_implicitdynamicsoccupancy.R 
##
## Purpose of script: Building implicit dynamics (multi-season) occupancy models in R
##
## Author: James E Paterson
##
## Date Created: 2024-06-02
##
## ---------------------------

## ----loaddata----------------------------------------------------

library(dplyr) # for data organization
library(ggplot2) # for plotting
# install.packages('RPresence', repo='https://www.mbr-pwrc.usgs.gov/mbrCRAN') # first time only
library(RPresence)

# Load detection history (100 sites with 15 visits each)
detection_history <- read.csv("dynamic_detection_history.csv", 
                              # First variable ("X") has row.names but not data
                              row.names = "X") 

# Examine data
head(detection_history)

# Load covariate data
effort <- read.csv("dynamic_effort.csv",
                   # First variable ("X") has row.names but not data
                   row.names = "X") 
site_cov <- read.csv("dynamic_site_cov.csv",
                     # First variable ("X") has row.names but not data
                     row.names = "X")

# RPresence uses createPao to create the input object for fitting occupancy models
frog_pao <- createPao( # data = a matrix with observed detection history 
                                          # 0's and 1's, one row per site, one column per survey
                                      data = as.matrix(detection_history),
                                      # nsurveyseason is a vector where the length = number of primary surveys
                                      # and the values are the number of secondary surveys
                                      nsurveyseason = rep(3,5),
                                      # survcov = survey covariates; a list where each element 
                                      # is a matrix the same size as "data"
                                      survcov = list(effort = as.matrix(effort)),
                                      # seasonncov is used for Explicit models (which we will compare to)
                                      seasncov = list(wetland_explicit_t = as.matrix(site_cov[,c(2, 5:8)])),
                                      # unitcov = dataframe with site rows x column variables, the site covariates
                                      # the first column of site_cov is the site number
                                      # Used for Explicit models (which we will compare to)
                                      unitcov = site_cov[,2:4],
                                      unitnames = as.character(site_cov$sites),
                                      title = "Frogs example") 

# Summary of pao occupancy model data
summary(frog_pao)

# Setting-up wetland change through time as an occupancy covariate for implicit dynamics occupancy model
# All site covariates used for Psi model must be included in this dataframe
# From occMod() documentation for implicit dynamics ("occMod_DO4"):
# psi.cov should be a data frame containing the unit-specific covariates to use for the occupancy
# component of the model, with number of rows = data$nunits or data$nunits*data$nseasons. If the 
# shorter version of the data frame is supplied the rows are recycled to the longer length.
implicit_occ_psi_covs <- data.frame(wetland_t = c(site_cov$wetland, site_cov$wetland.2,
                               site_cov$wetland.3, site_cov$wetland.4,
                               site_cov$wetland.5),
                               temp = rep(site_cov$temp, 5),
                               prec = rep(site_cov$prec, 5))

## ----buildoccMod-----------------------------------------------------------------------------------

# Fit models with occMod()
# type = "do.4" for implicit dynamics occupancy
implicit_occ_m1 <- occMod(model = list(psi ~wetland_t + temp + prec,
                                       p ~ effort), 
                          cov.list = list(psi.cov = implicit_occ_psi_covs),
                          data = frog_pao,
                          type = "do.4")

# Look at occupancy coefficients
implicit_occ_m1$beta$psi

# Look at detection coefficients
implicit_occ_m1$beta$p

## ----derivedparameters-----------------------------------------------------------------------------

# Show first 6 rows of derived colonization (gamma) parameter
# row names = site_primary survey period
head(implicit_occ_m1$derived$gamma)

# Show first 6 rows of derived extinction (epsilon) parameter
head(implicit_occ_m1$derived$epsilon)

# Test whether epsion = 1 - gamma for first row
implicit_occ_m1$derived$epsilon$est[1] == 1 - implicit_occ_m1$derived$gamma$est[1]

## ----predictedoccupancy----------------------------------------------------------------------------

# Mean predicted occupancy probabilities can be accessed using:
predicted_occupancy <- implicit_occ_m1$real$psi %>%
  mutate(site_time = row.names(.),
         site = rep(x = c(1:100), 5),
         time = rep(x = 1:5, each = 100))

predicted_occupancy_summary <- predicted_occupancy %>%
  group_by(time) %>%
  summarize(mean_occ = mean(est),
            lower_0.95 = mean(est) - 1.96*mean(se),
            upper_0.95 = mean(est) + 1.96*mean(se)) %>%
  mutate(type = "implicit dynamics")

## ----explicitdynamicmodels, echo = FALSE, message = FALSE, warning = FALSE-------------------------

# Fit explicit dynamic model in RPresence
explicit_occ_m1 <- occMod(model = list(psi ~wetland + temp + prec, # wetland + temp + prec 
                                       gamma ~wetland_explicit_t,
                                       epsilon ~wetland_explicit_t,
                                       p ~ effort),
                          data = frog_pao,
                          type = "do.1")

# Mean predicted occupancy probabilities from the explicit dynamics model to compare
explicit_predicted_occupancy <- explicit_occ_m1$real$psi %>%
  rbind(., 
        explicit_occ_m1$derived$psi) %>%
  mutate(site_time = row.names(.),
         site = rep(x = c(1:100), 5),
         time = rep(x = 1:5, each = 100))

# Summarize predicted occupancy by time across sites
explicit_predicted_occupancy_summary <- explicit_predicted_occupancy %>%
  group_by(time) %>%
  summarize(mean_occ = mean(est),
            lower_0.95 = mean(est) - 1.96*mean(se),
            upper_0.95 = mean(est) + 1.96*mean(se)) %>%
  mutate(type = "explicit dynamics")

# Bind rows from predictions of implicit and explicit dynamic occupancy models
predicted_occupancy_summary <- predicted_occupancy_summary  %>%
  rbind(explicit_predicted_occupancy_summary)

## ----plotderivedoccupancy----

# Position dodge the 
pd = position_dodge(width = 0.1)

ggplot(predicted_occupancy_summary, 
       aes(x = time, y = mean_occ, group = type, col = type)) +
  geom_errorbar(aes(ymin = lower_0.95,
                    ymax = upper_0.95),
                width = 0,
                position = pd) +
  geom_point(position = pd) +
  geom_line(position = pd) +
  scale_y_continuous(limits = c(0,1.0),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  scale_colour_manual(values = c("#7570B3", "#66A61E")) +
  labs(x = "Year", y = "Predicted occupancy probability") +
  theme_classic() +
  theme(legend.position = c(0.75, 0.75))
