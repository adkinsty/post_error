library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(shinystan)

setwd('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/memory_color_long/')
source('../../../forced_response/func.R')

dat <- read_csv("dat.csv")

results <- rpm(
  dat = dat, 
  mod_type = 'rpm1',
  re_run = TRUE, 
  R_mu = 1, 
  R_sigma = 1, 
  R_beta = 1,
  cov_mu = c("err0"), 
  cov_sigma = c("err0"), 
  cov_beta = c("err0"))

# follow up modeling of post-error data

dat_err0 <- dat %>% filter(err0 == .5)

results <- rpm(
  dat = dat_err0, 
  exp_name = 'memory_err0',
  mod_type = 'rpm1',
  re_run = TRUE)

results <- rpm(
  dat = dat_err0, 
  mod_type = 'rpm_mix',
  exp_name = 'memory_err0',
  re_run = TRUE, 
  R_mu = 2, 
  R_sigma = 2)

exp_name <- 'memory_err0'
fit1 <- read_rds(sprintf('fit/rpm1_%s_mu_1_sigma_1_beta_1_alpha_1.rds',exp_name))
fit_mix <- read_rds(sprintf('fit/rpm_mix_%s_mu_2_sigma_2_beta_1_alpha_1_lambda_1.rds',exp_name))

shinystan::launch_shinystan(fit_mix)

rnorm(1e4, -2, .5) %>%  inv_logit() %>% qplot() + xlim(0, 1)
rnorm(1e4, -1, .5) %>% inv_logit() %>% qplot() + xlim(0, 1)
rnorm(1e4, 0, .5) %>% inv_logit() %>% qplot() + xlim(0, 1)
rcauchy(1e4, 0, 5) %>%  inv_logit() %>% qplot() + xlim(0, 1)
rnorm(1e4, 1, .5) %>% inv_logit() %>% qplot() + xlim(0, 1)
