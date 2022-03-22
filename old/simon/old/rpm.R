library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686

setwd('/Users/adkinsty/Dropbox (University of Michigan)/from_box/side_projects/distraction/post_error/simon/')
source("func.R")

dat <- read_csv("dat.csv")

fit <- rpm(data = dat, re_run = TRUE, 
  R_mu = 2, R_sigma = 2, R_beta = 2,
  cov_mu = c("err0"), cov_sigma = c("err0"), cov_beta = c())