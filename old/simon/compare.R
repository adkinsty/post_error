library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(loo)

exp_name <- 'simon'

setwd(sprintf('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/%s/', exp_name))
source("func.R")



dat <- read_csv("dat.csv")


mu_sigma <- read_rds(sprintf('loo/loo_rpm_%s_mu_2err0_sigma_2err0_beta_2_alpha_1.rds',exp_name)) 
mu_beta <- read_rds(sprintf('loo/loo_rpm_%s_mu_2err0_sigma_2_beta_2err0_alpha_1.rds',exp_name)) 
mu_sigma_beta <- read_rds(sprintf('loo/loo_rpm_%s_mu_2err0_sigma_2err0_beta_2err0_alpha_1.rds',exp_name)) 


loo_compare(mu_sigma, mu_beta, mu_sigma_beta)
