library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(tidybayes) # Kay M (2020). tidybayes: Tidy Data and Geoms for Bayesian Models. doi: 10.5281/zenodo.1308151, R package version 2.3.1, http://mjskay.github.io/tidybayes/.
library(bayestestR)
library(bayesplot)
library(brms)

exp_name <- 'memory_color_short'

setwd(sprintf('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/%s/', exp_name))
source('../../../forced_response/func.R')

# Data

mod_name <- 'rpm1_mu_1err0_sigma_1err0_beta_1err0'

dat <- read_csv('dat.csv')

fit <- read_rds(sprintf('fit/%s_%s_alpha_1.rds',exp_name, mod_name)) 

draws <- as_draws_df(fit)

get_mu(draws, R_p = 1, R_mu = 1, covar = FALSE) %>% spread_draws(mu) %>% median_qi()
get_sigma(draws, R_p = 1, R_sigma = 1, covar = FALSE) %>% spread_draws(sigma) %>% median_qi()

get_mu(draws, R_p = 1, R_mu = 1, covar = TRUE) %>% spread_draws(mu[v]) %>% median_qi()
get_sigma(draws, R_p = 1, R_sigma = 1, covar = TRUE) %>% spread_draws(sigma[v]) %>% median_qi()

plt <- plot_prep1_cov(draws, exp_name)

delta_mu <- extract(fit, 'delta_mu[2]') %>% unlist() 
  delta_mu %>% p_direction() * sign(median(delta_mu))

delta_sigma <- extract(fit, 'delta_sigma[2]') %>% unlist() 
  delta_sigma %>% p_direction() * sign(median(delta_sigma))

delta_beta <- extract(fit, 'delta_beta[2]') %>% unlist() 
  delta_beta %>% p_direction() * sign(median(delta_beta))

plt_dat <- dat %>% smooth_obs(cov = c('err0')) %>% mutate(err0 = cov1)
plt_dat %>% 
  ggplot(aes(x = t, y = obs_mu, 
    ymax = obs_mu + obs_se * 1.96, ymin = obs_mu - obs_se * 1.96, 
    colour = as.character(err0), fill = as.character(err0), group = err0)) +
  sw_geom + my_theme + 
  scale_colour_manual(values = WONG[c(4, 7)]) +
  scale_fill_manual(values = WONG[c(4, 7)]) +
  ggsave(sprintf('fig/%s_obs.jpg', exp_name), units = 'in', height = 2, width = 3, dpi = 1000)

plt_dat <- dat %>% smooth_obs(cov = c('err0', 'err00')) %>% mutate(err0 = cov1, err00 = cov2)
plt_dat %>% 
  ggplot(aes(x = t, y = obs_mu, 
    ymax = obs_mu + obs_se * 1.96, ymin = obs_mu - obs_se * 1.96, 
    colour = as.character(err00), fill = as.character(err00), group = err00)) +
  sw_geom + my_theme +
  scale_colour_manual(values = WONG[c(4, 7)]) +
  scale_fill_manual(values = WONG[c(4, 7)]) +
  facet_wrap(err0 ~ ., labeller = labeller(err0 = c('-0.5' = 'After correct', '0.5' = 'After error'))) +
  ggsave(sprintf('fig/%s_err00_obs.jpg', exp_name), units = 'in', height = 2, width = 5.5, dpi = 1000)

plt_dat <- dat %>% smooth_obs_pred(fit = fit, cov = c('err0')) %>% mutate(err0 = cov1)
plt_dat %>% 
  ggplot(aes(x = t, y = pred, 
    ymax = upper, ymin = lower, colour = as.character(err0), fill = as.character(err0), group = err0)) +
  sw_geom + my_theme + 
  scale_y_continuous('Predicted Accuracy', breaks = seq(0, 1, .25), labels = c('0', '.25', '.5', '.75', '1')) + 
  scale_colour_manual(values = WONG[c(4, 7)]) +
  scale_fill_manual(values = WONG[c(4, 7)]) +
  ggsave(sprintf('fig/%s_pred.jpg', mod_name), units = 'in', height = 2, width = 3, dpi = 1000)

tmp <- plt_dat %>% 
  mutate(obs = obs_mu) %>%
  pivot_longer(cols = c(obs, pred), names_to = 'resp', values_to = 'mu') %>%
  mutate(upper = ifelse(resp == 'obs', obs_mu + obs_se * 1.96, upper),
    lower = ifelse(resp == 'obs', obs_mu - obs_se * 1.96, lower))
tmp %>%  
  ggplot(aes(x = t, y = mu, 
    ymax = upper, ymin = lower, alpha = as.character(resp), fill = as.character(resp), colour = as.character(resp))) +
  sw_geom +
  scale_colour_manual(values = WONG[c(1, 6)]) +
  scale_fill_manual(values = WONG[c(1, 6)]) +
  facet_wrap(err0 ~ ., labeller = labeller(err0 = c('-0.5' = 'After correct', '0.5' = 'After error'))) +
  scale_y_continuous('Accuracy', breaks = seq(0, 1, .25), labels = c('0', '.25', '.5', '.75', '1')) + 
  my_theme + ggsave(sprintf('fig/%s_obs_pred.jpg', mod_name), units = 'in', height = 2, width = 5.5, dpi = 1000)


color_scheme_set('purple')

hist_theme <- theme(axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.line.x = element_line(size = .25),
    axis.ticks.x = element_line(size = .25),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 8))

# intercepts
mcmc_hist(fit, pars = c('delta_mu[1]'), binwidth = .005, 
  transform = list('delta_mu[1]' = 'inv_logit')) + 
  my_theme + 
  geom_histogram(aes(x = inv_logit(rnorm(8e3, -1, .5))), alpha = .2, binwidth = .005) +
  scale_x_continuous(expression(mu), breaks = seq(0, .6, .2)) +
  coord_cartesian(xlim = c(0, .6), ylim = c(0, 4000)) + 
  hist_theme +
  ggsave(sprintf('fig/%s_mu_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(fit, pars = c('delta_sigma[1]'), binwidth = .005, 
  transform = list('delta_sigma[1]' = 'inv_logit')) + 
  my_theme + 
  geom_histogram(aes(x = inv_logit(rnorm(8e3, -2, .5))), alpha = .2, binwidth = .005) +
  scale_x_continuous(expression(sigma), breaks = c(0, .1, .2, .3, .4)) +
  coord_cartesian(xlim = c(0, .4), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_sigma_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(fit, pars = c('delta_beta[1]'), binwidth = .005, 
  transform = list('delta_beta[1]' = 'inv_logit')) + 
  my_theme + 
  geom_histogram(aes(x = inv_logit(rnorm(8e3, 1, .5))), alpha = .2, binwidth = .005) +
  scale_x_continuous(expression(beta), breaks = c(.6, .7, .8, .9, 1)) +
  coord_cartesian(xlim = c(.6, 1), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_beta_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

# reward differences
mu_diff <- get_diff(get_mu(draws, R_p = 1, R_mu = 1, covar = TRUE), par = 'mu', R = 1)
sigma_diff <- get_diff(get_sigma(draws, R_p = 1, R_sigma = 1, covar = TRUE), par = 'sigma', R = 1)
beta_diff <- get_diff(get_beta(draws, R_p = 1, R_beta = 1, covar = TRUE), par = 'beta', R = 1)

prior_start_mu <- rnorm(8e3, -1, .5)
prior_slope <- rnorm(8e3, 0, .5)

mcmc_hist(mu_diff, pars = c('mu_diff'), binwidth = .005) + 
  my_theme + 
  geom_histogram(aes(x = inv_logit(prior_start_mu + prior_slope * .5) - inv_logit(prior_start_mu + prior_slope * -.5)), alpha = .2, binwidth = .005) +
  scale_x_continuous(expression(Delta[mu]), breaks = seq(-.2, .2, .1)) +
  coord_cartesian(xlim = c(-.2, .2), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_diff_mu_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

prior_start_sig <- rnorm(8e3, -2, .5)

mcmc_hist(sigma_diff, pars = c('sigma_diff'), binwidth = .005) + 
  my_theme + 
  geom_histogram(aes(x = inv_logit(prior_start_sig + prior_slope * .5) - inv_logit(prior_start_sig + prior_slope * -.5)), alpha = .2, binwidth = .005) +
  scale_x_continuous(expression(Delta[sigma]), breaks = seq(-.2, .2, .1)) +
  coord_cartesian(xlim = c(-.2, .2), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_diff_sigma_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

prior_start_bet <- rnorm(8e3, 1, .5)

mcmc_hist(beta_diff, pars = c('beta_diff'), binwidth = .005) + 
  my_theme + 
  geom_histogram(aes(x = inv_logit(prior_start_bet + prior_slope * .5) - inv_logit(prior_start_bet + prior_slope * -.5)), alpha = .2, binwidth = .005) +
  scale_x_continuous(expression(Delta[beta]), breaks = seq(-.2, .2, .1)) +
  coord_cartesian(xlim = c(-.2, .2), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_diff_beta_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)














exp_name <- 'memory_err0'
fit1 <- read_rds(sprintf('fit/%s_rpm1_mu_1_sigma_1_beta_1_alpha_1.rds',exp_name))
fit_mix <- read_rds(sprintf('fit/%s_rpm_mix_mu_2_sigma_2_beta_1_alpha_1.rds',exp_name))

draws1 <- as_draws_df(fit1)
draws_mix <- as_draws_df(fit_mix)

plot_prep1(draws1, exp_name)

plot_prep_mix(draws_mix, exp_name)

get_lambda(draws_mix) %>% spread_draws(lambda) %>% median_qi()
get_mu(draws_mix, R_p = 2, R_mu = 2, covar = FALSE) %>% spread_draws(mu[r]) %>% median_qi()
get_sigma(draws_mix, R_p = 2, R_sigma = 2, covar = FALSE) %>% spread_draws(sigma[r]) %>% median_qi()

