library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(tidybayes) # Kay M (2020). tidybayes: Tidy Data and Geoms for Bayesian Models. doi: 10.5281/zenodo.1308151, R package version 2.3.1, http://mjskay.github.io/tidybayes/.
library(bayestestR)
library(bayesplot)
library(brms)

setwd('/Users/adkinsty/Dropbox (University of Michigan)/from_box/side_projects/distraction/post_error/simon/')
source("func.R")

# Data

exp_name <- 'arrow'
mod_name <- 'skew'

dat <- read_csv("dat.csv")

fit <- read_rds(sprintf('fit/rpm_%s_mu_2err0_sigma_2err0_beta_2_alpha_1.rds',exp_name)) 
fit <- read_rds(sprintf('fit/rpm_%s_%s_mu_2err0_sigma_2err0_lambda_2err0_beta_2_alpha_1.rds', mod_name, exp_name)) 
draws <- as_draws_df(fit)

get_mu(draws, covar = FALSE) %>% spread_draws(mu[r]) %>% median_qi()
get_sigma(draws, covar = FALSE) %>% spread_draws(sigma[r]) %>% median_qi()
get_lambda(draws, covar = FALSE) %>% spread_draws(lambda[r]) %>% median_qi()
get_beta(draws, covar = FALSE) %>% spread_draws(beta[r]) %>% median_qi()


get_mu(draws, covar = TRUE) %>% spread_draws(mu[r, v]) %>% median_qi()
get_sigma(draws, covar = TRUE) %>% spread_draws(sigma[r,v]) %>% median_qi()
get_lambda(draws, covar = TRUE) %>% spread_draws(lambda[r, v]) %>% median_qi()


#plt <- plot_prep(draws, exp_name)
plt <- plot_prep_skew(draws, exp_name)


plt <- plot_mu(draws, exp_name)
plt <- plot_beta(draws, exp_name)
plt <- plot_sigma(draws, exp_name)



extract(fit, 'delta_mu[2,1]') %>% unlist() %>% p_direction()
extract(fit, 'delta_mu[2,2]') %>% unlist() %>% p_direction()

extract(fit, 'delta_sigma[2,1]') %>% unlist() %>% p_direction()
extract(fit, 'delta_sigma[2,2]') %>% unlist() %>% p_direction()

extract(fit, 'delta_lambda[2,1]') %>% unlist() %>% p_direction()
extract(fit, 'delta_lambda[2,2]') %>% unlist() %>% p_direction()

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
    facet_wrap(stim_rep ~ .) +
  sw_geom + my_theme + 
  ggsave(sprintf('fig/%s_obs_stim_rep.jpg', exp_name), units = 'in', height = 2.5, width = 8, dpi = 1000)

plt_dat <- dat %>% smooth_obs_pred(fit = fit)
plt_dat %>% 
  ggplot(aes(x = t, y = pred, 
    ymax = upper, ymin = lower,
    colour = factor(congr), fill = factor(congr), alpha = as.character(err0), group = interaction(congr, err0))) +
  sw_geom + scale_y_continuous("Predicted Accuracy", breaks = seq(0, 1, .25), labels = c("0", ".25", ".5", ".75", "1")) + 
  my_theme + ggsave(sprintf('fig/%s_%s_pred.jpg', exp_name, mod_name), units = 'in', height = 2, width = 3, dpi = 1000)

mcmc_scatter(fit, pars = c("delta_mu[2,1]", "delta_mu[2,2]"), alpha = .2, size = .2) + 
  my_theme + 
  scale_x_continuous(expression(Delta[mu[h]]), breaks = c(-1, 0, 1)) +
  scale_y_continuous(expression(Delta[mu[g]]), breaks = c(-1, 0, 1)) +
  geom_smooth(method = "lm",se = FALSE, size = .25, colour = "black") +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  ggsave(sprintf('fig/%s_%s_mcmc_scatter.jpg', exp_name, mod_name), units = 'in', height = 2, width = 2.1, dpi = 1000)
  
mcmc_scatter(fit, pars = c("delta_sigma[2,1]", "delta_sigma[2,2]"), alpha = .1, size = .1) + 
  my_theme + 
  scale_x_continuous(expression(Delta[sigma[h]]), breaks = c(-2, 0, 2)) +
  scale_y_continuous(expression(Delta[sigma[g]]), breaks = c(-2, 0, 2)) +
  geom_smooth(method = "lm",se = FALSE, size = .25, colour = "black") +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  ggsave(sprintf('fig/%s_%s_delta_sigma_mcmc_scatter.jpg', exp_name, mod_name), units = 'in', height = 2, width = 2, dpi = 1000)

mcmc_scatter(fit, pars = c("delta_lambda[2,1]", "delta_lambda[2,2]"), alpha = .1, size = .1) + 
  my_theme + 
  scale_x_continuous(expression(Delta[lambda[h]]), breaks = c(-2, 0, 2)) +
  scale_y_continuous(expression(Delta[lambda[g]]), breaks = c(-2, 0, 2)) +
  geom_smooth(method = "lm",se = FALSE, size = .25, colour = "black") +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  ggsave(sprintf('fig/%s_%s_delta_lambda_mcmc_scatter.jpg', exp_name, mod_name), units = 'in', height = 2, width = 2, dpi = 1000)


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
