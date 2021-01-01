## ---------------------------
##
## Script name: 04_dynamicoccupancy.R 
##
## Purpose of script: Building dynamic (multi-season) occupancy models in R
##
## Author: Dr. James Paterson
##
## Date Created: 2021-01-01
##
## 
## Email: james.earle.paterson@gmail.com
##
## ---------------------------


## ----loaddata-------------------------------------------
library(dplyr) # for data organization
library(magrittr) # for piping
library(ggplot2) # for plotting
# install.packages("unmarked") # first time only
library(unmarked) # for occupancy models

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

# Build an unmarkedmultFramOccu
frog_unmarkedMultFrame <- unmarkedMultFrame( # y is a matrix with observed detection history 
                                          # 0's and 1's, one row per site, one column per survey
                                      y = as.matrix(detection_history),
                                      # numPrimary is the number of primary surveys; 
                                      # function derives # of secondary surveys based on
                                      # ncol(y) / numPrimary
                                      numPrimary = 5,
                                      # obsCovs = observation covariates in a list, 
                                      # each variable has site rows x survey columns
                                      obsCovs = list(effort = effort),
                                      # siteCovs = dataframe with site rows x column variables
                                      # the first column of site_cov is the site number
                                      # the 5th:8th columns are used below in yearlySiteCovs
                                      siteCovs = site_cov[,2:4],
                                      # "yearlySiteCovs" are a list for site variables that differ between primary survey periods. 
                                      # In our example, the amount of wetland changes between primary survey periods and 
                                      # This data is stored within site_cov. There must be one column per primary survey (5 in our example)
                                      yearlySiteCovs = list(wetland_time = site_cov[,c(2, 5:8)])) 

# S4 class for colext occupancy model data
summary(frog_unmarkedMultFrame)

## ----buildcolext---------------------------------------------------------
dynamic_occ_m1 <- colext(
  # Psi depends on initial habitat, climate
  psiformula= ~wetland + temp + prec,
  # colonization depends on wetland change
  gammaformula = ~wetland_time, 
  # extinction depends on wetland change
  epsilonformula = ~wetland_time,
  # detection depends on survey effort
  pformula = ~effort,
  # data must be a unmarkedMultFrame
  data = frog_unmarkedMultFrame,
  # method is optim method, leave as "BFGS"
  method = "BFGS")

## ----dynamicgof----------
# Mackenzie-Bailey GOF test
# Simulate capture history data (if model correct). Compare obs X2 to sim X2
# Must simulate 1000-5000 times for good distribution estimates
# Likely to take a Very Long Time on large data sets
# with nsim = 10, takes my Mackbook ~ 60 seconds
mb.boot <- AICcmodavg::mb.gof.test(dynamic_occ_m1,
                                   nsim = 10) # Must be much higher than five to be useful
print(mb.boot, digit.vals = 4, digits.chisq = 4)


## ----modelsummary--------------------------------------------------------
summary(dynamic_occ_m1)

## ----sitepredictions-------------------
# Get the smoothed occupancy probabilities for site 25 
# [2, , 25] = 2 for occupied, " " for all 5 surveys, 25 for site 25; remove [] to print whole array
dynamic_occ_m1@smoothed[2, , 25]

## ----smoothedmeanpredictions---------------------------------------------
# Mean smoothed occupancy probabilities can be accessed using:
# dynamic_occ_m1@smoothed.mean, or
# smoothed(dynamic_occ_m1)

# Calculate SE for derived occupancy predictions using bootstrap, larger B requires longer time
m1 <- nonparboot(dynamic_occ_m1, 
                 B = 10)

# Predicted occupancy in each year (with SE)
predicted_occupancy <- data.frame(year = c(1:5),
                             smoothed_occ = smoothed(dynamic_occ_m1)[2,],
                             SE = m1@smoothed.mean.bsse[2,])

## ----plotderivedoccupancy----
ggplot(predicted_occupancy, 
       aes(x = year, y = smoothed_occ)) +
  geom_errorbar(aes(ymin = smoothed_occ-SE,
                    ymax = smoothed_occ+SE),
                width = 0) +
  geom_point(size = 3) +
  geom_line() +
  scale_y_continuous(limits = c(0,1.0),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = "Year", y = "Smoothed occupancy (derived)") +
  theme_classic()