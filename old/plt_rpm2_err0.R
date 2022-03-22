library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(tidybayes) # Kay M (2020). tidybayes: Tidy Data and Geoms for Bayesian Models. doi: 10.5281/zenodo.1308151, R package version 2.3.1, http://mjskay.github.io/tidybayes/.
library(bayestestR)
library(bayesplot)
library(brms)


proj_name <- 'post_error'
exp_name <- 'flanker'

setwd(sprintf('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/%s/%s/', proj_name, exp_name))
source('../../../forced_response/func.R')

# Data

mod_name <- 'rpm2_hier_mu_err0_sigma_err0_beta_err0'

dat <- read_csv('dat.csv')

fit <- read_rds(sprintf('fit/%s_%s.rds',exp_name, mod_name)) 
fit_loo <- read_rds(sprintf('loo/loo_%s_%s.rds',exp_name, mod_name))
draws <- as_draws_df(fit)

get_mu(draws, R_mu = 2, covar = FALSE) %>% spread_draws(mu[r]) %>% median_qi()
get_sigma(draws, R_sigma = 2, covar = FALSE) %>% spread_draws(sigma[r]) %>% median_qi()
get_beta(draws, R_beta = 2, covar = FALSE) %>% spread_draws(beta[r]) %>% median_qi()

get_mu(draws, R_mu = 2, covar = TRUE) %>% spread_draws(mu[r,v]) %>% median_qi()
get_sigma(draws, R_sigma = 2, covar = TRUE) %>% spread_draws(sigma[r,v]) %>% median_qi()
get_beta(draws, R_beta = 2, covar = TRUE) %>% spread_draws(beta[r,v]) %>% median_qi()

plt <- plot_prep2_cov(draws, mod_name)

delta_mu_h <- extract(fit, 'delta_mu[2,1]') %>% unlist() 
  delta_mu_h %>% p_direction() * sign(median(delta_mu_h))
delta_mu_g <- extract(fit, 'delta_mu[2,2]') %>% unlist() 
  delta_mu_g %>% p_direction() * sign(median(delta_mu_g))

delta_sigma_h <- extract(fit, 'delta_sigma[2,1]') %>% unlist() 
  delta_sigma_h %>% p_direction() * sign(median(delta_sigma_h))
delta_sigma_g <- extract(fit, 'delta_sigma[2,2]') %>% unlist() 
  delta_sigma_g %>% p_direction() * sign(median(delta_sigma_g))

delta_beta_h <- extract(fit, 'delta_beta[2,1]') %>% unlist() 
delta_beta_h %>% p_direction() * sign(median(delta_beta_h))
delta_beta_g <- extract(fit, 'delta_beta[2,2]') %>% unlist() 
delta_beta_g %>% p_direction() * sign(median(delta_beta_g))

# dpd <- sw_ttest(tmp, cov = 'congr')
# dpd_c <- sw_ttest(tmp %>% filter(congr > 0))
# dpd_i <- sw_ttest(tmp %>% filter(congr < 0))

# dpd_dat <- tibble(time = 50:700, pd = abs(dpd_c), congr = 'c') %>% 
#   bind_rows(tibble(time = 50:700, pd = abs(dpd_i), congr = 'i'))
# write_csv(dpd_dat, sprintf('all/tab/%s_dpd.csv', exp_name))

# dpd_dat %>%
#   mutate(congr = factor(congr, c('i', 'c'))) %>%
#   ggplot(aes(x = time, y = pd, colour = congr, alpha = pd > .9)) + 
#   pd_geom + my_theme + 
#   ggsave(sprintf('all/fig/%s_pd.jpg', exp_name), units = 'in', height = 1.5, width = 4, dpi = 1000)

plt_dat <- dat %>% 
  mutate(X = PT, Y = accuracy) %>%
  smooth_obs(covar = c('congr','err0'), width = 100) %>% 
  mutate(congr = covar1, err0 = covar2)

plt_dat %>% 
  ggplot(aes(x = PT, y = obs_mu, 
    ymax = obs_mu + obs_se * 1.96, ymin = obs_mu - obs_se * 1.96, 
    colour = interaction(congr, err0), fill = interaction(congr, err0), group = interaction(congr, err0))) +
  sw_geom + my_theme + #raw + 
  ggsave('fig/obs.jpg', units = 'in', height = 2, width = 3, dpi = 1000)

plt_dat <- dat %>%
  mutate(X = PT, Y = accuracy) %>%
  smooth_obs_pred(fit = fit, covar = c('congr','err0'), width = 100) %>% 
  mutate(congr = covar1, err0 = covar2)
plt_dat %>% 
  ggplot(aes(x = PT, y = pred, 
    ymax = upper, ymin = lower, 
    colour = interaction(congr, err0), fill = interaction(congr, err0), group = interaction(congr, err0))) +
  sw_geom + scale_y_continuous('Predicted Accuracy', breaks = seq(0, 1, .25), labels = c('0', '.25', '.5', '.75', '1')) + 
  my_theme + ggsave(sprintf('fig/%s_pred.jpg', mod_name), units = 'in', height = 2, width = 3, dpi = 1000)

tmp <- plt_dat %>% 
  mutate(obs = obs_mu) %>%
  pivot_longer(cols = c(obs, pred), names_to = 'resp', values_to = 'mu') %>%
  mutate(upper = ifelse(resp == 'obs', obs_mu + obs_se * 1.96, upper),
    lower = ifelse(resp == 'obs', obs_mu - obs_se * 1.96, lower))
tmp %>%  
  ggplot(aes(x = PT, y = mu, 
    ymax = upper, ymin = lower, alpha = as.character(resp), fill = as.character(resp), colour = as.character(resp))) +
  sw_geom +
  scale_fill_manual(values = WONG[c(1, 3)]) +
  scale_colour_manual(values = WONG[c(1, 3)]) +
  facet_wrap(congr ~ err0, 
    labeller = labeller(congr = c('-0.5' = 'Incongruent', '0.5' = 'Congruent'), err0 = c('-0.5' = 'Post correct', '0.5' = 'Post error'))) +
  scale_y_continuous('Accuracy', breaks = seq(0, 1, .25), labels = c('0', '.25', '.5', '.75', '1')) + 
  my_theme + ggsave(sprintf('fig/%s_obs_pred.jpg', mod_name), units = 'in', height = 4, width = 6, dpi = 1000)


color_scheme_set('purple')

hist_theme <- theme(axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.line.x = element_line(size = .25),
    axis.ticks.x = element_line(size = .25),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 8))

# intercepts
mcmc_hist(fit, pars = c('delta_mu[1,1]'), binwidth = .005, 
  transform = list('delta_mu[1,1]' = 'inv_logit')) + 
  my_theme + 
  scale_x_continuous(expression(mu^h), breaks = c(0, .1, .2, .3, .4)) +
  coord_cartesian(xlim = c(0, .4), ylim = c(0, 4000)) + 
  hist_theme +
  ggsave(sprintf('fig/%s_mu_h_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(fit, pars = c('delta_mu[1,2]'), binwidth = .005, 
  transform = list('delta_mu[1,2]' = 'inv_logit')) + 
  my_theme + 
  scale_x_continuous(expression(mu^g), breaks = c(0, .1, .2, .3, .4)) +
  coord_cartesian(xlim = c(0, .4), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_mu_g_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(fit, pars = c('delta_sigma[1,1]'), binwidth = .005, 
  transform = list('delta_sigma[1,1]' = 'inv_logit')) + 
  my_theme + 
  scale_x_continuous(expression(sigma^h), breaks = c(0, .1, .2, .3, .4)) +
  coord_cartesian(xlim = c(0, .4), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_sigma_h_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(fit, pars = c('delta_sigma[1,2]'), binwidth = .005, 
  transform = list('delta_sigma[1,2]' = 'inv_logit')) + 
  my_theme + 
  scale_x_continuous(expression(sigma^g), breaks = c(0, .1, .2, .3, .4)) +
  coord_cartesian(xlim = c(0, .4), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_sigma_g_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(fit, pars = c('delta_beta[1,1]'), binwidth = .005, 
  transform = list('delta_beta[1,1]' = 'inv_logit')) + 
  my_theme + 
  scale_x_continuous(expression(beta^h), breaks = c(.6, .7, .8, .9, 1)) +
  coord_cartesian(xlim = c(.6, 1), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_beta_h_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(fit, pars = c('delta_beta[1,2]'), binwidth = .005, 
  transform = list('delta_beta[1,2]' = 'inv_logit')) + 
  my_theme + 
  scale_x_continuous(expression(beta^g), breaks = c(.6, .7, .8, .9, 1)) +
  coord_cartesian(xlim = c(.6, 1), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_beta_g_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

# err0 differences
mu_diff <- get_diff(get_mu(draws, covar = TRUE), par = 'mu')
sigma_diff <- get_diff(get_sigma(draws, covar = TRUE), par = 'sigma')
beta_diff <- get_diff(get_beta(draws, covar = TRUE), par = 'beta')

mcmc_hist(mu_diff, pars = c('mu_diff[h]'), binwidth = .005) + 
  my_theme + 
  scale_x_continuous(expression(Delta[mu]^h), breaks = c(-.1, 0, .1)) +
  coord_cartesian(xlim = c(-.1, .1), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_diff_mu_h_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(mu_diff, pars = c('mu_diff[g]'), binwidth = .005) + 
  my_theme + 
  scale_x_continuous(expression(Delta[mu]^g), breaks = c(-.1, 0, .1)) +
  coord_cartesian(xlim = c(-.1, .1), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_diff_mu_g_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(sigma_diff, pars = c('sigma_diff[h]'), binwidth = .005) + 
  my_theme + 
  scale_x_continuous(expression(Delta[sigma]^h), breaks = c(-.1, 0, .1)) +
  coord_cartesian(xlim = c(-.1, .1), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_diff_sigma_h_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(sigma_diff, pars = c('sigma_diff[g]'), binwidth = .005) + 
  my_theme + 
  scale_x_continuous(expression(Delta[sigma]^g), breaks = c(-.1, 0, .1)) +
  coord_cartesian(xlim = c(-.1, .1), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_diff_sigma_g_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(beta_diff, pars = c('beta_diff[h]'), binwidth = .005) + 
  my_theme + 
  scale_x_continuous(expression(Delta[beta]^h), breaks = c(-.1, 0, .1)) +
  coord_cartesian(xlim = c(-.1, .1), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_diff_beta_h_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

mcmc_hist(beta_diff, pars = c('beta_diff[g]'), binwidth = .005) + 
  my_theme + 
  scale_x_continuous(expression(Delta[beta]^g), breaks = c(-.1, 0, .1)) +
  coord_cartesian(xlim = c(-.1, .1), ylim = c(0, 4000)) +
  hist_theme + 
  ggsave(sprintf('fig/%s_diff_beta_g_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)



dat %>% 
  mutate(congr = factor(congr), err0 = factor(err0)) %>%
  ggplot(aes(x = PT, fill = interaction(congr, err0), colour = interaction(congr, err0), group = interaction(congr, err0))) +
  geom_histogram(position = "identity", bins = 100, alpha = .5) +
  scale_x_continuous('PT (ms)', breaks = seq(0, 1000, 100)) +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 100)) +
  scale_colour_manual(values = c(WONG[2], WONG[3], WONG[7], WONG[6])) +
  scale_fill_manual(values = c(WONG[2], WONG[3], WONG[7], WONG[6])) +
  my_theme +
  hist_theme +
  ggsave(sprintf('fig/%s_PT_hist_err0.jpg', exp_name), units = 'in', height = 2, width = 3, dpi = 1000)
