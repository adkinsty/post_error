library(brms)
library(tidyverse)

setwd('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/simon/')
source('func.R')


true_dat <- read_csv('dat.csv')
all_dat <- read_csv('dat/all_test_data.csv')

all_dat %>% 
  filter(PT > 0 & PT < 1000) %>%
  ggplot(aes(x = PT, group = trial_resp.corr, colour = trial_resp.corr, fill = trial_resp.corr)) +
  geom_density(alpha = .5) + 
  facet_wrap(accuracy ~.)

exp_name <- 'fake'

#skews <- rep(seq(-4, 5, 1), each = 2)
#skews <- rep(c(0, 2), each = 2)

N <- 1e4 #nrow(true_dat)

dat <- tibble( 
    exp = exp_name,
    congr = rep(c(-.5, .5), N / 2), #true_dat$congr,
    err0 = rep(c(-.5, -.5, .5, .5), N / 4), #true_dat$err0, 
    time = runif(N), # true_dat$time, 
    PT =  time*1000, #true_dat$PT,
    lambda = ifelse(err0 > 0, 10, 0), #rep(skews, N / 4),
    alpha = .5,
    beta_h = ifelse(congr > 0, .85, 1 - .85), 
    beta_g = .95,
    mu_h = .20,
    mu_g = .35,
    sigma_h = .05,
    sigma_g = .10,
    lambda_h = lambda,
    lambda_g = lambda,
    phi_h = pskew_normal(q = time, mu = mu_h, sigma = sigma_h, alpha = lambda_h), 
    phi_g = pskew_normal(q = time, mu = mu_g, sigma = sigma_g, alpha = lambda_g), 
    p_guess = (1 - phi_h) * (1 - phi_g) * alpha,
    p_habit = phi_h * (1 - phi_g) * beta_h, 
    p_goal = (1 - phi_h) * phi_g * beta_g, 
    p_both = phi_h * phi_g * beta_g, 
    theta = p_guess + p_habit + p_goal + p_both) %>%
  rowwise() %>%
  mutate(accuracy = rbernoulli(n = 1, p = theta))


plt_dat <- dat %>% group_by(congr, lambda) %>% summarise(time = 50:700, obs = sw_smooth(x = PT, y = accuracy), pred = sw_smooth(x = PT, y = theta))

plt_dat %>% 
  ggplot(aes(x = time, y = obs, ymax = 0, ymin = 0, colour = factor(congr), fill = factor(congr), alpha = factor(lambda))) +
  sw_geom + scale_y_continuous("Observed Accuracy", breaks = seq(0, 1, .25), labels = c("0", ".25", ".5", ".75", "1")) + 
  #scale_alpha_manual(values = seq(.25, .75, .05)) + geom_line(data = plt_dat %>% filter(lambda == 0), aes(x = t, y = pred), colour = 'black', alpha = 1, size = .25) +
  my_theme + ggsave('fig/skew_sim_ideal_obs.jpg', units = 'in', height = 2, width = 3, dpi = 1000)
plt_dat %>% 
  ggplot(aes(x = time, y = pred, ymax = 0, ymin = 0, colour = factor(congr), fill = factor(congr), alpha = factor(lambda))) +
  sw_geom + scale_y_continuous("Predicted Accuracy", breaks = seq(0, 1, .25), labels = c("0", ".25", ".5", ".75", "1")) + 
  #scale_alpha_manual(values = seq(.25, .75, .05)) + geom_line(data = plt_dat %>% filter(lambda == 0), aes(x = t, y = pred), colour = 'black', alpha = 1, size = .25) +
  my_theme + ggsave('fig/skew_sim_ideal_pred.jpg', units = 'in', height = 2, width = 3, dpi = 1000)

fit <- rpm_skew(data = dat, re_run = TRUE, 
  R_mu = 2, R_sigma = 2, R_lambda = 2, R_beta = 2,
  cov_mu = c(), cov_sigma = c(), 
  cov_lambda = c("err0"), cov_beta = c())

#fit <- read_rds(sprintf('fit/rpm_skew_%s_mu_2_sigma_2_lambda_2err0_beta_2_alpha_1.rds',exp_name)) 
fit <- read_rds(sprintf('fit/rpm_skew_%s_mu_2err0_sigma_2err0_lambda_2err0_beta_2_alpha_1.rds',exp_name)) 

draws <- as_draws_df(fit)
get_mu(draws) %>% spread_draws(mu[r, v]) %>% median_qi()
get_sigma(draws) %>% spread_draws(sigma[r, v]) %>% median_qi()
get_lambda(draws) %>% spread_draws(lambda[r, v]) %>% median_qi()

get_mu(draws, covar = FALSE) %>% spread_draws(mu[r]) %>% median_qi()
get_sigma(draws, covar = FALSE) %>% spread_draws(sigma[r]) %>% median_qi()
get_lambda(draws, covar = FALSE) %>% spread_draws(lambda[r]) %>% median_qi()

get_beta(draws, covar = FALSE) %>% spread_draws(beta[r]) %>% median_qi()

draws %>% spread_draws(delta_lambda[r, v]) %>% median_qi()
draws %>% spread_draws(delta_sigma[r, v]) %>% median_qi()
draws %>% spread_draws(delta_mu[r, v]) %>% median_qi()
