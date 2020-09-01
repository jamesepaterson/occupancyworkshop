## ---------------------------
##
## Script name: 02_singleseasonoccupancy.R 
##
## Purpose of script: Introduction to occupancy models in R (single season)
##
## Author: Dr. James Paterson
##
## Date Created: 2020-09-01
##
## 
## Email: james.earle.paterson@gmail.com
##
## ---------------------------


## ----loaddata-------------------------------------------
# install.packages("unmarked") # first time only
library(unmarked)

# Load detection history (100 sites with 10 visits each)
detection_history <- read.csv("detection_history.csv", 
                              # First variable ("X") has row.names but not data
                              row.names = "X") 

# Examine data
head(detection_history)

# Create unmarkedFrameOccu that holds the data
sample.unmarkedFrame_simple <- unmarkedFrameOccu( # y is a matrix with observed detection history 
                                          # (0's and 1's, one row per site, one column per survey)
                                      y = as.matrix(detection_history)) 

# S4 class for occupancy model data
summary(sample.unmarkedFrame_simple)

## ----buildbasicoccu------------------------------------------------------
# Build basic single-season occupancy model with intercepts only (one estimate for detection, one for occupancy)
occu.m1 <- occu(formula = ~1 # detection formula first
                          ~1, # occupancy formula second, 
                  data = sample.unmarkedFrame_simple)

summary(occu.m1) # Show AIC, estimates (on logit scale), SE, z-scores

# To get real estimate of occupancy (with 95% CI)
predict(occu.m1, 
        newdata = data.frame(site = 1),
        type = "state")

# To get real estimate of detection (with 95% CI)
predict(occu.m1, 
        newdata = data.frame(site = 1),
        type = "det")

# Equivalent to inverse logit
boot::inv.logit(coef(occu.m1)[1]) # Real estimate of occupancy
boot::inv.logit(coef(occu.m1)[2]) # Real estimate of detection

## ----covariatesload------------------------------------------------------
# Load covariate data
effort <- read.csv("effort.csv",
                   # First variable ("X") has row.names but not data
                   row.names = "X") 
observers <- read.csv("observers.csv",
                      # First variable ("X") has row.names but not data
                      row.names = "X") 
site_cov <- read.csv("site_cov.csv",
                     # First variable ("X") has row.names but not data
                     row.names = "X") 

# Build a new unmarkedFramOccu
sample.unmarkedFrame_cov <- unmarkedFrameOccu( # y is a matrix with observed detection history 
                                          # (0's and 1's, one row per site, one column per survey)
                                      y = as.matrix(detection_history),
                                      # obsCovs = observation covariates in a list, 
                                      # each variable has site rows x survey columns
                                      obsCovs = list(effort = effort,
                                                     observers = observers),
                                      # siteCovs = dataframe with site rows x column variables
                                      siteCovs = site_cov) 

# S4 class for occupancy model data
summary(sample.unmarkedFrame_cov)

## ----buildoccucov--------------------------------------------------------
occu.m2 <- occu(formula = ~effort + observers # detection formula first
                          ~forest + agri, # occupancy formula second,
                data = sample.unmarkedFrame_cov)

# Summarize
summary(occu.m2)

## ----covpredict----------
# Predict effect on new data set to see how occupancy changes with `forest`
predict_m2_forest <- cbind(predict(occu.m2,
                             newdata = data.frame(forest = seq(min(site_cov$forest, 
                                                           na.rm = TRUE),
                                                       max(site_cov$forest, 
                                                           na.rm = TRUE), 
                                                       by = 0.01),
                                                  agri = mean(site_cov$agri)),
                             type = "state"),
                           data.frame(forest = seq(min(site_cov$forest, 
                                                           na.rm = TRUE),
                                                       max(site_cov$forest, 
                                                           na.rm = TRUE), 
                                                       by = 0.01),
                                                  agri = mean(site_cov$agri)))

# Predict effect on new data set to see how occupancy changes with `agri`
predict_m2_agri <- cbind(predict(occu.m2,
                             newdata = data.frame(agri = seq(min(site_cov$agri, 
                                                           na.rm = TRUE),
                                                       max(site_cov$agri, 
                                                           na.rm = TRUE), 
                                                       by = 0.01),
                                                  forest = mean(site_cov$forest)),
                             type = "state"),
                           data.frame(agri = seq(min(site_cov$agri, 
                                                           na.rm = TRUE),
                                                       max(site_cov$agri, 
                                                           na.rm = TRUE), 
                                                       by = 0.01),
                                                  forest = mean(site_cov$forest)))

## ----plotrelationships----
library(ggplot2)
# Plot relationship with forest
ggplot(data = predict_m2_forest, aes(x = forest, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray") +
  stat_smooth(method = "loess", col = "black", se = FALSE) +
  labs(x = "Forest (scaled)", y = "Predicted Occupancy Probability") +
  theme_classic()

# Plot relationship with agri
ggplot(data = predict_m2_agri, aes(x = agri, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray") +
  stat_smooth(method = "loess", col = "black", se = FALSE) +
  labs(x = "Agriculture (scaled)", y = "Predicted Occupancy Probability") +
  theme_classic()

## ----mbgof---------------------------------------------------------------
# install.packages("AICcmodavg") # First time only
library(AICcmodavg)

# Do Mackenzie-Bailey goodness of fit test for single-season occupancy model
m2_mb.gof.boot <- mb.gof.test(occu.m2,
                                # Demonstrate with small number of sims (10), 
                                # but then change to large number (e.g. 1000)
                                nsim = 50)

# View Results
m2_mb.gof.boot