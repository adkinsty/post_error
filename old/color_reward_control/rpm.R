library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(shinystan)

setwd('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/color_reward_control/')
source("func.R")

dat <- read_csv("dat.csv")


fit <- rpm_skew(data = dat, re_run = TRUE, 
  R_mu = 2, R_sigma = 2, R_lambda = 2, R_beta = 2,
  cov_mu = c("err0"), cov_sigma = c("err0"), 
  cov_lambda = c("err0"), cov_beta = c())


fit <- read_rds('fit/rpm_skew_color_reward_control_mu_2err0_sigma_2err0_lambda_2err0_beta_2_alpha_1.rds')
launch_shinystan(fit)

draws <- as_draws_df(fit)
get_mu(draws) %>% spread_draws(mu[r, v]) %>% median_qi()
get_sigma(draws) %>% spread_draws(sigma[r, v]) %>% median_qi()
get_lambda(draws) %>% spread_draws(lambda[r, v]) %>% median_qi()
get_lambda(draws, covar = FALSE) %>% spread_draws(lambda[r]) %>% median_qi()

draws %>% spread_draws(delta_lambda[k, r]) %>% median_qi()
draws %>% spread_draws(delta_sigma[k, r]) %>% median_qi()
draws %>% spread_draws(delta_mu[k, r]) %>% median_qi()
