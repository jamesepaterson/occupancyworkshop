## ---------------------------
##
## Script name: 03_singleseasonoccupancy_part2_modelcomparison.R 
##
## Purpose of script: Building and comparing occupancy models in R (single season)
##
## Author: Dr. James Paterson
##
## Date Created: 2020-11-08
##
## 
## Email: james.earle.paterson@gmail.com
##
## ---------------------------


## ----loaddata-------------------------------------------
# install.packages("unmarked") # first time only
library(unmarked)
library(dplyr)
library(magrittr)
library(ggplot2)

# Load detection history (100 sites with 10 visits each)
detection_history <- read.csv("detection_history.csv", 
                              # First variable ("X") has row.names but not data
                              row.names = "X") 

# Examine data
head(detection_history)

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

# Add a random variable to site_cov (to demonstrate that not all variables should be included in the model)
# We'll assume this variable measures the amount of wetland in some buffer around the site
# Because it is a random number, it shouldn't affect occupancy in our example.
site_cov$wetland <- rnorm(n = nrow(site_cov), mean = 0, sd = 4)

# Build an unmarkedFramOccu
sample.unmarkedFrame_cov <- unmarkedFrameOccu( # y is a matrix with observed detection history 
                                          # 0's and 1's, one row per site, one column per survey
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
# Fit general model all variables
occu_p_full_psi_full <- occu(formula = ~effort + observers # detection formula first
                          ~forest + agri + wetland, # occupancy formula second,
                data = sample.unmarkedFrame_cov)

## ----fitsubmodels--------
# Workflow:
# 1. fit general model
# 2. Fit all sub-sets of the general model
# 3. Make subset based on models that are well-supported (< 10 delta AICc from most supported)
# 4. Base inferenece on model-averaged predictions from supported models

# 1. General model fit above (occu_p_full_psi_full)

# 2. Fit all sub-sets of the general model
# install.packages("MuMIn") # first time only
library(MuMIn)
occu_dredge <- dredge(global.model = occu_p_full_psi_full,
                        # the rank argument can use AICc, QAICc, or others (see help page)
                        rank = "AICc")

# Look at model selection table (just showing top 5 models, only the top 2 have any support)
occu_dredge[1:5,]

# 3. Create model list for those that are relatively well-supported (e.g. < 10 delta AICc from most supported)
# In this case, there is only one additional model (compared to general model) with any support
# We will fit this model, then put the 2 supported models in a list
occu_p_full_psi_agr_frs <- occu(formula = ~effort + observers # detection formula first
                          ~forest + agri, # occupancy formula second,
                data = sample.unmarkedFrame_cov)

# Make model list
occu_model_list <- list(occ_full = occu_p_full_psi_full,
                        occ_sub  = occu_p_full_psi_agr_frs)


## ----modelaveragedpredictions----------
# With lots of models, it will be slow
library(AICcmodavg)
occu_modavg_psi_predict <- modavgPred(occu_model_list, 
                                        # c.hat = 1, # to change variance inflation factor, default = 1) 
                                        parm.type = "psi", # psi = occupancy, can also be "detect" for detection probability
                                        newdata = sample.unmarkedFrame_cov@siteCovs)[c("mod.avg.pred",
                                                                                       "lower.CL",
                                                                                       "upper.CL")]


## Put predictions, CI, and all site covariates into one data frame
occu_modavg_psi_predict_df <- data.frame(Predicted = occu_modavg_psi_predict$mod.avg.pred,
                                    lower = occu_modavg_psi_predict$lower.CL,
                                    upper = occu_modavg_psi_predict$upper.CL,
                                    site_cov)

# Look at first values
head(occu_modavg_psi_predict_df)

## ----mapoccu---------------------------------------
occu_modavg_psi_predict_df <- occu_modavg_psi_predict_df %>%
  mutate(x = rep(1:10, 10), # Adding pseudo-coordinates as if 100 sites are in 10 x 10 grid 
         y = rep(1:10, each = 10))

# Note, you will probably plot on map and use polygons, not points
# This is an example just to show how you can easily visualize the predicted occupancy across all sites
ggplot(data = occu_modavg_psi_predict_df, aes(x = x, y = y, fill = Predicted)) +
  geom_point(size = 10, pch = 22) +
  theme_classic() +
  scale_fill_gradient(low = "blue", high = "yellow", name = "Predicted\noccupancy") +
  ggtitle("Predicted occupancy")

## ----modelaveraging_covariates---------------------
# First, set-up a new dataframe to predict along a sequence of the covariate.
# Predicting requires all covariates, so let's hold the other covariates constant at their mean value
occu_forest_newdata <- data.frame(forest = seq(min(site_cov$forest), 
                                                 max(site_cov$forest), by = 0.5),
                                  agri = mean(site_cov$agri), # hold other variables constant
                                  wetland = mean(site_cov$wetland)) # hold other variables constant
  
# Model-averaged prediction of occupancy and confidence interval
occu_forest_pred <- modavgPred(occu_model_list,
                               # c.hat =    # to change variance inflation factor, default = 1) 
                               parm.type = "psi", # psi = occupancy
                               newdata = occu_forest_newdata)[c("mod.avg.pred",
                                                                "lower.CL",
                                                                "upper.CL")]

# Put prediction, confidence interval, and covariate values together in a data frame
occu_forest_pred_df <- data.frame(Predicted = occu_forest_pred$mod.avg.pred,
                                  lower = occu_forest_pred$lower.CL,
                                  upper = occu_forest_pred$upper.CL,
                                  occu_forest_newdata)

# Plot the relationship
occu_forest_pred_plot <- ggplot(occu_forest_pred_df, aes(x = forest, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, linetype = "dashed") +
  geom_path(size = 1) +
  labs(x = "Forest cover (standardized)", y = "Occupancy probability") +
  theme_classic() +
  coord_cartesian(ylim = c(0,1)) +
  theme(text = element_text(family = "HelveticaNeue", colour = "black"),
        axis.text = element_text(colour = "black"))
occu_forest_pred_plot

