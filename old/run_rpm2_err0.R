library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(shinystan)

proj_name <- 'post_error'
exp_name <- 'memory'
mod_name <- 'rpm2_hier'

setwd(sprintf('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/%s/%s/', proj_name, exp_name))
source('../../forced_response/func.R')

dat <- read_csv('dat.csv')

results <- fit_rpm(
  dat = dat, 
  exp_name = exp_name,
  mod_name = mod_name,
  prior_only = 0, 
  gq = 1,
  re_run = TRUE, 
  covar = 'err0')

# library(tidybayes)

# fit <- read_rds(sprintf('fit/%s_%s_mu_err0_sigma_err0.rds', exp_name, mod_name))
# draws <- as_draws_df(fit)








# posterior intervals

# mu_loc <- draws %>% spread_draws(mu_loc[r]) %>% mutate(mu_loc = inv_logit(mu_loc)) %>% median_qi()
# mu_loc %>% ggplot(aes(y = mu_loc, ymax = .upper, ymin = .lower, x = factor(r), colour = factor(r))) + 
#   geom_point(position = position_dodge(.3)) + 
#   geom_errorbar(position = position_dodge(.3), width = 0) + 
#   coord_flip(ylim = c(0, 1)) + 
#   my_theme + 
#   theme(
#     panel.grid.major.x = element_line(), 
#     panel.grid.minor.x = element_blank()) + 
#   ggsave('fig/tmp.jpg', dpi = 1000, height = 1, width = 3)

# sigma_loc <- draws %>% spread_draws(sigma_loc[r]) %>% mutate(sigma_loc = inv_logit(sigma_loc)) %>% median_qi()
# sigma_loc %>% ggplot(aes(y = sigma_loc, ymax = .upper, ymin = .lower, x = factor(r), colour = factor(r))) + 
#   geom_point(position = position_dodge(.3)) + 
#   geom_errorbar(position = position_dodge(.3), width = 0) + 
#   coord_flip(ylim = c(0, 1)) + 
#   my_theme + 
#   theme(
#     panel.grid.major.x = element_line(), 
#     panel.grid.minor.x = element_blank()) + 
#   ggsave('fig/tmp.jpg', dpi = 1000, height = 1, width = 3)


# beta_loc <- draws %>% spread_draws(beta_loc[r]) %>% mutate(beta_loc = inv_logit(beta_loc)) %>% median_qi()
# beta_loc %>% ggplot(aes(y = beta_loc, ymax = .upper, ymin = .lower, x = factor(r), colour = factor(r))) + 
#   geom_point(position = position_dodge(.3)) + 
#   geom_errorbar(position = position_dodge(.3), width = 0) + 
#   coord_flip(ylim = c(0, 1)) + 
#   my_theme + 
#   theme(
#     panel.grid.major.x = element_line(), 
#     panel.grid.minor.x = element_blank()) + 
#   ggsave('fig/tmp.jpg', dpi = 1000, height = 1, width = 3)

# alpha_loc <- draws %>% spread_draws(alpha_loc) %>% mutate(alpha_loc = inv_logit(alpha_loc)) %>% median_qi()
# alpha_loc %>% ggplot(aes(y = alpha_loc, ymax = .upper, ymin = .lower, x = 1)) + 
#   geom_point(position = position_dodge(.3)) + 
#   geom_errorbar(position = position_dodge(.3), width = 0) + 
#   coord_flip(ylim = c(0, 1)) + 
#   my_theme + 
#   theme(
#     panel.grid.major.x = element_line(), 
#     panel.grid.minor.x = element_blank()) + 
#   ggsave('fig/tmp.jpg', dpi = 1000, height = 1, width = 3)

# mu0 <- draws %>% spread_draws(mu0[j, r]) %>% mutate(mu0 = inv_logit(mu0)) %>% median_qi()
# mu0 %>% ggplot(aes(y = mu0, ymax = .upper, ymin = .lower, x = factor(j), colour = factor(r))) + 
#   geom_point(position = position_dodge(.3)) + 
#   geom_errorbar(position = position_dodge(.3), width = 0) + 
#   coord_flip(ylim = c(0, 1)) + 
#   my_theme + 
#   theme(
#     panel.grid.major.x = element_line(), 
#     panel.grid.minor.x = element_blank()) + 
#   ggsave('fig/tmp.jpg', dpi = 1000, height = 3, width = 3)

# sigma0 <- draws %>% spread_draws(sigma0[j, r]) %>% mutate(sigma0 = inv_logit(sigma0)) %>% median_qi()
# sigma0 %>% ggplot(aes(y = sigma0, ymax = .upper, ymin = .lower, x = factor(j), colour = factor(r))) + 
#   geom_point(position = position_dodge(.3)) + 
#   geom_errorbar(position = position_dodge(.3), width = 0) + 
#   coord_flip(ylim = c(0, 1)) + 
#   my_theme +
#   theme(
#     panel.grid.major.x = element_line(), 
#     panel.grid.minor.x = element_blank()) +  
#   ggsave('fig/tmp.jpg', dpi = 1000, height = 3, width = 3)

# beta0 <- draws %>% spread_draws(beta0[j, r]) %>% mutate(beta0 = inv_logit(beta0)) %>% median_qi()
# beta0 %>% ggplot(aes(y = beta0, ymax = .upper, ymin = .lower, x = factor(j), colour = factor(r))) + 
#   geom_point(position = position_dodge(.3)) + 
#   geom_errorbar(position = position_dodge(.3), width = 0) + 
#   coord_flip(ylim = c(0, 1)) + 
#   my_theme + 
#   theme(
#     panel.grid.major.x = element_line(), 
#     panel.grid.minor.x = element_blank()) + 
#   ggsave('fig/tmp.jpg', dpi = 1000, height = 3, width = 3)

# alpha0 <- draws %>% spread_draws(alpha0[j]) %>% mutate(alpha0 = inv_logit(alpha0)) %>% median_qi()
# alpha0 %>% ggplot(aes(y = alpha0, ymax = .upper, ymin = .lower, x = factor(j))) + 
#   geom_point(position = position_dodge(.3)) + 
#   geom_errorbar(position = position_dodge(.3), width = 0) + 
#   coord_flip(ylim = c(0, 1)) + 
#   my_theme + 
#   theme(
#     panel.grid.major.x = element_line(), 
#     panel.grid.minor.x = element_blank()) + 
#   ggsave('fig/tmp.jpg', dpi = 1000, height = 3, width = 3)

# delta_mu <- draws %>% spread_draws(delta_mu[j, r]) %>% median_qi()
# delta_mu %>% ggplot(aes(y = delta_mu, ymax = .upper, ymin = .lower, x = factor(j), colour = factor(r))) + 
#   geom_point(position = position_dodge(.3)) + 
#   geom_errorbar(position = position_dodge(.3), width = 0) + 
#   coord_flip(ylim = c(-2, 2)) + 
#   my_theme + 
#   theme(
#     panel.grid.major.x = element_line(), 
#     panel.grid.minor.x = element_blank()) + 
#   ggsave('fig/tmp.jpg', dpi = 1000, height = 3, width = 3)

# delta_sigma <- draws %>% spread_draws(delta_sigma[j, r]) %>% median_qi()
# delta_sigma %>% ggplot(aes(y = delta_sigma, ymax = .upper, ymin = .lower, x = factor(j), colour = factor(r))) + 
#   geom_point(position = position_dodge(.3)) + 
#   geom_errorbar(position = position_dodge(.3), width = 0) + 
#   coord_flip(ylim = c(-2, 2)) + 
#   my_theme + 
#   theme(
#     panel.grid.major.x = element_line(), 
#     panel.grid.minor.x = element_blank()) + 
#   ggsave('fig/tmp.jpg', dpi = 1000, height = 3, width = 3)
