library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(tidybayes) # Kay M (2020). tidybayes: Tidy Data and Geoms for Bayesian Models. doi: 10.5281/zenodo.1308151, R package version 2.3.1, http://mjskay.github.io/tidybayes/.
library(bayestestR)
library(bayesplot)
library(brms)

exp_name <- 'simon'

setwd(sprintf('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/%s/', exp_name))
source("func.R")

# Data

mod_name <- 'mu_2err0_sigma_2err0_beta_2err0'

dat <- read_csv("dat.csv")

fit <- read_rds(sprintf('fit/rpm_%s_%s_alpha_1.rds',exp_name, mod_name)) 

draws <- as_draws_df(fit)

get_mu(draws, covar = FALSE) %>% spread_draws(mu[r]) %>% median_qi()
get_sigma(draws, covar = FALSE) %>% spread_draws(sigma[r]) %>% median_qi()
#get_lambda(draws, covar = FALSE) %>% spread_draws(lambda[r]) %>% median_qi()
get_beta(draws, covar = FALSE) %>% spread_draws(beta[r]) %>% median_qi()


get_mu(draws, covar = TRUE) %>% spread_draws(mu[r, v]) %>% median_qi()
get_sigma(draws, covar = TRUE) %>% spread_draws(sigma[r,v]) %>% median_qi()
#get_lambda(draws, covar = TRUE) %>% spread_draws(lambda[r, v]) %>% median_qi()


plt <- plot_prep(draws, exp_name)
#plt <- plot_prep_skew(draws, exp_name)


plt <- plot_mu(draws, exp_name)
plt <- plot_beta(draws, exp_name)
plt <- plot_sigma(draws, exp_name)



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

#extract(fit, 'delta_lambda[2,1]') %>% unlist() %>% p_direction()
#extract(fit, 'delta_lambda[2,2]') %>% unlist() %>% p_direction()

deltas <- delta_summary(fit, exp_name)

mu_diff1 <- get_mu_diff(get_mu(draws))
mu_diff1$`mu_diff[hi]` %>% quantile(probs = c(.05, .5, .95))
mu_diff1$`mu_diff[lo]` %>% quantile(probs = c(.05, .5, .95))

mu_diff2 <- get_mu_diff2(mu_diff1)
mu_diff2$`mu_diff2` %>% quantile(probs = c(.05, .5, .95))
mu_diff2$`mu_diff2` %>% p_direction()



plt_dat <- dat %>% smooth_obs()
plt_dat %>% 
  ggplot(aes(x = t, y = obs_mu, 
    ymax = obs_mu + obs_se, ymin = obs_mu - obs_se,
    colour = factor(congr), fill = factor(congr), alpha = as.character(err0), group = interaction(congr, err0))) +
    #facet_wrap(stim_rep ~ .) +
  sw_geom + my_theme + 
  ggsave(sprintf('fig/%s_obs.jpg', exp_name), units = 'in', height = 2, width = 3, dpi = 1000)

plt_dat <- dat %>% smooth_obs_pred(fit = fit)
plt_dat %>% 
  ggplot(aes(x = t, y = pred, 
    ymax = upper, ymin = lower,
    colour = factor(congr), fill = factor(congr), alpha = as.character(err0), group = interaction(congr, err0))) +
  sw_geom + scale_y_continuous("Predicted Accuracy", breaks = seq(0, 1, .25), labels = c("0", ".25", ".5", ".75", "1")) + 
  my_theme + ggsave(sprintf('fig/%s_pred.jpg', mod_name), units = 'in', height = 2, width = 3, dpi = 1000)


plt_dat <- dat %>% smooth_obs_pred(fit = fit)
tmp <- plt_dat %>% 
  mutate(obs = obs_mu) %>%
  pivot_longer(cols = c(obs, pred), names_to = 'resp', values_to = 'mu') %>%
  mutate(upper = ifelse(resp == 'obs', obs_mu + obs_se, upper),
    lower = ifelse(resp == 'obs', obs_mu - obs_se, lower))
tmp %>%  
  ggplot(aes(x = t, y = mu, 
    ymax = upper, ymin = lower,
    colour = factor(congr), fill = factor(congr), alpha = as.character(resp), group = interaction(congr, resp, err0))) +
  sw_geom +
  facet_wrap(err0 ~ ., labeller = labeller(err0 = c('-0.5' = 'After correct', '0.5' = 'After error'))) +
  scale_y_continuous("Accuracy", breaks = seq(0, 1, .25), labels = c("0", ".25", ".5", ".75", "1")) + 
  my_theme + ggsave(sprintf('fig/%s_obs_pred.jpg', mod_name), units = 'in', height = 2, width = 5, dpi = 1000)


mcmc_scatter(fit, pars = c("delta_mu[2,1]", "delta_mu[2,2]"), alpha = .2, size = .2) + 
  my_theme + 
  scale_x_continuous(expression(Delta[mu[h]]), breaks = c(-1, 0, 1)) +
  scale_y_continuous(expression(Delta[mu[g]]), breaks = c(-1, 0, 1)) +
  geom_smooth(method = "lm",se = FALSE, size = .25, colour = "black") +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  ggsave(sprintf('fig/%s_delta_mu_mcmc_scatter.jpg', mod_name), units = 'in', height = 2, width = 2.1, dpi = 1000)
  
mcmc_scatter(fit, pars = c("delta_sigma[2,1]", "delta_sigma[2,2]"), alpha = .1, size = .1) + 
  my_theme + 
  scale_x_continuous(expression(Delta[sigma[h]]), breaks = c(-2, 0, 2)) +
  scale_y_continuous(expression(Delta[sigma[g]]), breaks = c(-2, 0, 2)) +
  geom_smooth(method = "lm",se = FALSE, size = .25, colour = "black") +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  ggsave(sprintf('fig/%s_delta_sigma_mcmc_scatter.jpg', mod_name), units = 'in', height = 2, width = 2.1, dpi = 1000)

mcmc_scatter(fit, pars = c("delta_beta[2,1]", "delta_beta[2,2]"), alpha = .1, size = .1) + 
  my_theme + 
  scale_x_continuous(expression(Delta[beta[h]]), breaks = c(-2, 0, 2)) +
  scale_y_continuous(expression(Delta[beta[g]]), breaks = c(-2, 0, 2)) +
  geom_smooth(method = "lm",se = FALSE, size = .25, colour = "black") +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  ggsave(sprintf('fig/%s_delta_beta_mcmc_scatter.jpg', mod_name), units = 'in', height = 2, width = 2.1, dpi = 1000)


# mcmc_scatter(fit, pars = c("delta_lambda[2,1]", "delta_lambda[2,2]"), alpha = .1, size = .1) + 
#   my_theme + 
#   scale_x_continuous(expression(Delta[lambda[h]]), breaks = c(-2, 0, 2)) +
#   scale_y_continuous(expression(Delta[lambda[g]]), breaks = c(-2, 0, 2)) +
#   geom_smooth(method = "lm",se = FALSE, size = .25, colour = "black") +
#   coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
#   ggsave(sprintf('fig/%s_delta_lambda_mcmc_scatter.jpg', mod_name), units = 'in', height = 2, width = 2.1, dpi = 1000)


dat %>%
  mutate(congr = as.factor(congr), err0 = as.factor(err0)) %>%
  group_by(congr, err0) %>% 
  summarise(
    mean = unlist(mean_cl_boot(y, na.rm = T)['y']), 
    lower = unlist(mean_cl_boot(y, na.rm = T)['ymin']), 
    upper = unlist(mean_cl_boot(y, na.rm = T)['ymax'])) %>%
  ggplot(
    aes(x = err0, y = mean, ymin = lower, ymax = upper, fill = congr, colour = congr)) +
  geom_col(position = position_dodge(.75), width = .725, alpha = .5) + 
  geom_errorbar(position = position_dodge(.75), width = .25) +
  scale_y_continuous('Accuracy', limits = c(0, 1)) +
  scale_x_discrete("Previous Trial Error", labels = c("0", "1")) +
  scale_fill_manual("Congruency", labels = c("i", "c"), values = WONG[c(7, 6)]) +
  scale_colour_manual("Congruency", labels = c("i", "c"), values = WONG[c(7, 6)]) +
  my_theme + 
  theme(panel.grid.major.x = element_blank()) +
  ggsave(sprintf('fig/%s_obs_col.jpg', exp_name), units = 'in', height = 2.25, width = 2.25, dpi = 1000)
