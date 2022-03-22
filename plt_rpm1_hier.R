library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(tidybayes) # Kay M (2020). tidybayes: Tidy Data and Geoms for Bayesian Models. doi: 10.5281/zenodo.1308151, R package version 2.3.1, http://mjskay.github.io/tidybayes/.
library(bayestestR)
library(bayesplot)
library(brms)
library(bayestestR)

proj_name <- 'post_error'
exp_name <- 'memory_color_long'
covar_name <- 'err0'
mod_type <- 'rpm1_hier'

setwd(sprintf('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/%s/%s/', proj_name, exp_name))
source('../../../forced_response/func.R')

# Data
dat <- read_csv('dat.csv') %>% mutate(PT_true = PT * 2)



# performance regression models
my_priors <- c(
  set_prior('normal(0, 1)', class = 'b'), 
  set_prior('normal(0, 1)', class = 'Intercept'), 
  set_prior('normal(0, 1)', class = 'sd'))
m <- read_rds('fit/brm.rds')
m <- brm(
  data = dat,
  formula = y ~ err0 + (err0 | participant), 
  family = 'bernoulli', 
  prior = my_priors,
  chains = 4, cores = 1, iter = 4000, warmup = 2000)
write_rds(m, 'fit/brm.rds')

m_early <- read_rds('fit/brm_early.rds')
m_early <- update(m, newdata = dat %>% filter(PT_true < 500), cores = 1)
write_rds(m_early, 'fit/brm_early.rds')

m_later <- read_rds('fit/brm_later.rds')
m_later <- update(m, newdata = dat %>% filter(PT_true > 1000), cores = 1)
write_rds(m_later, 'fit/brm_later.rds')





# response preparation model
mod_name <- sprintf('%s_mu_%s_sigma_%s_beta_%s', mod_type, covar_name, covar_name, covar_name)
fit <- read_rds(sprintf('fit/%s_%s.rds',exp_name, mod_name)) 
draws <- as_draws_df(fit)

p <- plot_prep1_cov(draws, exp_name, transform = inv_logit2, time = seq(0, 1.5, .001)) +  
  ggsave(filename = sprintf('fig/%s_prep_dens1_cov.jpg', mod_name), units = 'in', height = 1.5, width = 3.5, dpi = 1000)

# p <- plot_prep1(draws, exp_name, transform = inv_logit2, time = seq(0, 1.5, .001)) +  
#   scale_x_continuous('Time (ms)', breaks = seq(0, 1500, 150), labels = as.character(seq(0, 1500, 150))) +
#   coord_cartesian(xlim = c(0, 1500), y = c(0, 10)) +   
#   ggsave(filename = sprintf('fig/%s_prep_dens1.jpg', mod_name), units = 'in', height = 1.5, width = 3.5, dpi = 1000)

get_mu(draws, R_mu = 1, covar = FALSE, transform = inv_logit2) %>% spread_draws(mu) %>% median_qi()
get_sigma(draws, R_sigma = 1, covar = FALSE, transform = inv_logit2) %>% spread_draws(sigma) %>% median_qi()
get_beta(draws, R_beta = 1, covar = FALSE, transform = inv_logit) %>% spread_draws(beta) %>% median_qi()
get_beta(draws, R_beta = 1, covar = TRUE) %>% spread_draws(beta[v]) %>% median_qi()


# group intercepts
mu_prior <- 2*inv_logit(rnorm(8e3, -.5, .5))

h <- mcmc_hist(fit, 
  pars = c('mu_loc'), 
  binwidth = .01, 
  transform = list('mu_loc' = 'inv_logit2'))
p <- h + 
  my_theme + 
  hist_theme +
  scale_x_continuous(expression(mu), breaks = seq(0, 1, .1)) +
  coord_cartesian(xlim = c(0, .5), ylim = c(0, 4000)) + 
  ggsave(sprintf('fig/%s_mu_mcmc_hist.jpg', mod_name), units = 'in', height = 2, width = 3, dpi = 1000)

h <- mcmc_hist(fit, pars = c('sigma_loc'), binwidth = .01, 
  transform = list('sigma_loc' = 'inv_logit2')) 
p <- h + 
  my_theme + 
  hist_theme + 
  geom_histogram(aes(x = 2*inv_logit(rnorm(8e3, -2, .5))), alpha = .4, binwidth = .005) +
  scale_x_continuous(expression(sigma), breaks = seq(0, 2, .2)) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 4000)) +
  ggsave(sprintf('fig/%s_sigma_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

p <- mcmc_hist(fit, pars = c('beta_loc'), binwidth = .005, 
  transform = list('beta_loc' = 'inv_logit')) + 
  my_theme + 
  hist_theme + 
  geom_histogram(aes(x = inv_logit(rnorm(8e3, 2, .5))), alpha = .4, binwidth = .005) +
  scale_x_continuous(expression(beta), breaks = c(.6, .7, .8, .9, 1)) +
  coord_cartesian(xlim = c(.6, 1), ylim = c(0, 4000)) +
  ggsave(sprintf('fig/%s_beta_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)


# group reward differences

mus <- get_mu(draws, R_mu = 1, covar = TRUE) 
sigmas <- get_sigma(draws, R_sigma = 1, covar = TRUE)
betas <- get_beta(draws, R_beta = 1, covar = TRUE)

mu_diff <- get_diff(mus, par = 'mu', R = 1) %>% subset_draws(variable = 'mu_diff')
mu_diff %>% median_qi()
mu_diff %>% p_direction()

sigma_diff <- get_diff(sigmas, par = 'sigma', R = 1) %>% subset_draws(variable = 'sigma_diff')
sigma_diff %>% median_qi()
sigma_diff %>% p_direction()

beta_diff <- get_diff(betas, par = 'beta', R = 1) %>% subset_draws(variable = 'beta_diff')
beta_diff %>% median_qi()
beta_diff %>% p_direction()

prior_start_mu <- rnorm(8e3, -.5, .5)
prior_slope <- rnorm(8e3, 0, .5)

x2 <- function(x) {
  return(2*x)
}

p <- mcmc_hist(mu_diff, pars = c('mu_diff'), binwidth = .005) + 
  my_theme + 
  hist_theme + 
  # geom_histogram(aes(x = 2*inv_logit(prior_start_mu + prior_slope * .5) - 2*inv_logit(prior_start_mu + prior_slope * -.5)), alpha = .4, binwidth = .005) +
  scale_x_continuous(expression(Delta[mu]), breaks = seq(-.2, .2, .1)) +
  coord_cartesian(xlim = c(-.15, .15), ylim = c(0, 4000)) +
  ggsave(sprintf('fig/%s_diff_mu_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1, dpi = 1000)

prior_start_sig <- rnorm(8e3, -2, .5)

p <- mcmc_hist(sigma_diff, pars = c('sigma_diff'), binwidth = .005) + 
  my_theme + 
  hist_theme + 
  # geom_histogram(aes(x = 2*inv_logit(prior_start_sig + prior_slope * .5) - 2*inv_logit(prior_start_sig + prior_slope * -.5)), alpha = .4, binwidth = .005) +
  scale_x_continuous(expression(Delta[sigma]), breaks = seq(-.2, .2, .1)) +
  coord_cartesian(xlim = c(-.15, .15), ylim = c(0, 4000)) +
  ggsave(sprintf('fig/%s_diff_sigma_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1, dpi = 1000)

prior_start_bet <- rnorm(8e3, 2, .5)

p <- mcmc_hist(beta_diff, pars = c('beta_diff'), binwidth = .005) + 
  my_theme + 
  hist_theme + 
  #geom_histogram(aes(x = inv_logit(prior_start_bet + prior_slope * .5) - inv_logit(prior_start_bet + prior_slope * -.5)), alpha = .4, binwidth = .005) +
  scale_x_continuous(expression(Delta[beta]), breaks = seq(-.2, .2, .1)) +
  coord_cartesian(xlim = c(-.15, .15), ylim = c(0, 4000)) +
  ggsave(sprintf('fig/%s_diff_beta_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1, dpi = 1000)
 

p <- mcmc_hist(beta_diff, pars = c('beta_diff'), binwidth = .001) + 
  my_theme + 
  hist_theme + 
  #geom_histogram(aes(x = inv_logit(prior_start_bet + prior_slope * .5) - inv_logit(prior_start_bet + prior_slope * -.5)), alpha = .4, binwidth = .005) +
  scale_x_continuous(expression(Delta[beta]), breaks = seq(-.04, .04, .02), labels = c('-.04', '-.02', '0', '.02', '.04')) +
  coord_cartesian(xlim = c(-.03, .03), ylim = c(0, 4000)) +
  ggsave(sprintf('fig/%s_diff_beta_mcmc_hist.jpg', mod_name), units = 'in', height = 1, width = 1, dpi = 1000)


# conditional accuracy
plt_dat <- dat %>% 
  mutate(X = PT, Y = accuracy) %>%
  smooth_obs_pred(fit, covar_names = c(covar_name))
# plt_dat <- dat %>% 
#   mutate(X = PT, Y = accuracy) %>%
#   smooth_obs(covar_names = c(covar_name), width = 100)
p <- plt_dat %>% 
  ggplot(aes(x = PT, y = obs_mu, ymax = obs_mu + obs_se * 1.96, ymin = obs_mu - obs_se * 1.96, colour = as.character(covar1), fill = as.character(covar1), group = covar1)) +
  sw_geom + 
  my_theme +  
  scale_x_continuous("PT (ms)", breaks = seq(0, 1000, 100), labels = as.character(seq(0, 2000, 200))) +
  scale_colour_manual(values = WONG[c(4, 7)]) +
  scale_fill_manual(values = WONG[c(4, 7)]) +
  ggsave(sprintf('fig/%s_obs.jpg', exp_name), units = 'in', height = 2, width = 3.5, dpi = 1000)

p <- plt_dat %>% 
  ggplot(aes(x = PT, y = pred, ymax = upper, ymin = lower, colour = as.character(covar1), fill = as.character(covar1), group = covar1)) +
  sw_geom + 
  my_theme + 
  scale_x_continuous("PT (ms)", breaks = seq(0, 1000, 100), labels = as.character(seq(0, 2000, 200))) +
  scale_y_continuous('Predicted Accuracy', breaks = seq(0, 1, .25), labels = c('0', '.25', '.5', '.75', '1')) + 
  scale_colour_manual(values = WONG[c(4, 7)]) +
  scale_fill_manual(values = WONG[c(4, 7)]) +
  ggsave(sprintf('fig/%s_pred.jpg', mod_name), units = 'in', height = 2, width = 3.5, dpi = 1000)

tmp <- plt_dat %>% 
  mutate(obs = obs_mu) %>%
  pivot_longer(cols = c(obs, pred), names_to = 'resp', values_to = 'mu') %>%
  mutate(upper = ifelse(resp == 'obs', obs_mu + obs_se * 1.96, upper),
    lower = ifelse(resp == 'obs', obs_mu - obs_se * 1.96, lower))
  
facet_labels <- if (covar_name == 'rew') {c('$1', '$2')} else { if (covar_name == 'congr0') {c('After incongruent', 'After congruent')} else {c('After correct', 'After error')}}

p <- tmp %>%  
  ggplot(aes(x = PT, y = mu, ymax = upper, ymin = lower, alpha = as.character(resp), fill = as.character(resp), colour = as.character(resp))) +
  sw_geom +
  scale_colour_manual(values = WONG[c(1, 3)]) +
  scale_fill_manual(values = WONG[c(1, 3)]) +
  scale_x_continuous("PT (ms)", breaks = seq(0, 1000, 100), labels = as.character(seq(0, 2000, 200))) +
  facet_wrap(covar1 ~ ., labeller = labeller(covar1 = c('-0.5' = facet_labels[1], '0.5' = facet_labels[2]))) +
  scale_y_continuous('Accuracy', breaks = seq(0, 1, .25), labels = c('0', '.25', '.5', '.75', '1')) + 
  my_theme + ggsave(sprintf('fig/%s_obs_pred.jpg', mod_name), units = 'in', height = 2, width = 6.5, dpi = 1000)

plt_dat <- dat %>% 
  mutate(X = PT, Y = accuracy) %>%
  smooth_obs(covar_names = c(covar_name, 'key_stim_rep'), width = 100)
p <- plt_dat %>% 
  ggplot(aes(x = PT, y = obs_mu, 
    ymax = obs_mu + obs_se * 1.96, ymin = obs_mu - obs_se * 1.96, 
    colour = as.character(covar1), fill = as.character(covar1), group = covar1)) +
  sw_geom + my_theme +  
  scale_x_continuous("PT (ms)", breaks = seq(0, 1000, 100), labels = as.character(seq(0, 2000, 200))) +
  scale_colour_manual(values = WONG[c(4, 7)]) +
  scale_fill_manual(values = WONG[c(4, 7)]) +
  facet_wrap(covar2 ~ ., labeller = labeller(covar2 = c('TRUE' = 'Goal is previous key', 'FALSE' = 'Goal is not previous key'))) +
  ggsave(sprintf('fig/%s_obs_rep.jpg', exp_name), units = 'in', height = 2, width = 6.5, dpi = 1000)


plt_dat <- draws %>% subset_draws('mu0') %>% 
  spread_draws(mu0[j]) %>% 
  mutate(mu0 = inv_logit2(mu0) * 1000) %>%
  median_qi() %>%
  ungroup() %>%
  arrange(mu0) %>%
  mutate(j_asc = factor(1:n()))
plt_dat %>%
  ggplot(aes(x = mu0, y = j_asc, xmax = .upper, xmin = .lower)) + 
  geom_pointinterval(size = .1) + 
  my_theme + 
  scale_x_continuous(expression(mu[0]^j)) +
  scale_y_discrete('Participant (j)', labels = rep('', nrow(plt_dat))) +
  coord_cartesian(xlim = c(0, 1000)) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()) + 
  ggsave(sprintf('fig/%s_mu0_intervals.jpg', mod_name), units = 'in', height = 3, width = 2, dpi = 1000)








# subject level intervals
plt_dat <- draws %>% subset_draws('mu0') %>% 
  spread_draws(mu0[j]) %>% 
  mutate(mu0 = inv_logit2(mu0) * 1000) %>%
  median_qi() %>%
  ungroup() %>%
  arrange(mu0) %>%
  mutate(j_asc = factor(1:n()))
plt_dat %>%
  ggplot(aes(x = mu0, y = j_asc, xmax = .upper, xmin = .lower)) + 
  geom_pointinterval(size = .1) + 
  my_theme + 
  scale_x_continuous(expression(mu[0]^j)) +
  scale_y_discrete('Participant (j)', labels = rep('', nrow(plt_dat))) +
  coord_cartesian(xlim = c(0, 1000)) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()) + 
  ggsave(sprintf('fig/%s_mu0_intervals.jpg', mod_name), units = 'in', height = 4, width = 2, dpi = 1000)

plt_dat <- draws %>% subset_draws('delta_mu') %>% 
  spread_draws(delta_mu[j]) %>% 
  median_qi() %>%
  ungroup() %>%
  arrange(delta_mu) %>%
  mutate(j_asc = factor(1:n()))
plt_dat %>%
  ggplot(aes(x = delta_mu, y = j_asc, xmax = .upper, xmin = .lower)) + 
  geom_pointinterval(size = .1) + 
  my_theme + 
  scale_x_continuous(expression(Delta[mu]^j)) +
  scale_y_discrete('Participant (j)', labels = rep('', nrow(plt_dat))) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()) + 
  ggsave(sprintf('fig/%s_delta_mu_intervals.jpg', mod_name), units = 'in', height = 4, width = 2, dpi = 1000)

plt_dat <- draws %>% subset_draws('beta0') %>% 
  spread_draws(beta0[j]) %>% 
  mutate(beta0 = inv_logit(beta0)) %>%
  median_qi() %>%
  ungroup() %>%
  arrange(beta0) %>%
  mutate(j_asc = factor(1:n()))
plt_dat %>%
  ggplot(aes(x = beta0, y = j_asc, xmax = .upper, xmin = .lower)) + 
  geom_pointinterval(size = .1) + 
  my_theme + 
  scale_x_continuous(expression(beta[0]^j)) +
  scale_y_discrete('Participant (j)', labels = rep('', nrow(plt_dat))) +
  coord_cartesian(xlim = c(0, 1)) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()) + 
  ggsave(sprintf('fig/%s_beta0_intervals.jpg', mod_name), units = 'in', height = 4, width = 2, dpi = 1000)

plt_dat <- draws %>% subset_draws('delta_beta') %>% 
  spread_draws(delta_beta[j]) %>% 
  median_qi() %>%
  ungroup() %>%
  arrange(delta_beta) %>%
  mutate(j_asc = factor(1:n()))
plt_dat %>%
  ggplot(aes(x = delta_beta, y = j_asc, xmax = .upper, xmin = .lower)) + 
  geom_pointinterval(size = .1) + 
  my_theme + 
  scale_x_continuous(expression(Delta[mu]^j)) +
  scale_y_discrete('Participant (j)', labels = rep('', nrow(plt_dat))) +
  #coord_cartesian(xlim = c(0, 1000)) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()) + 
  ggsave(sprintf('fig/%s_delta_beta_intervals.jpg', mod_name), units = 'in', height = 4, width = 2, dpi = 1000)


# proportion response type
# tmp_dat <- dat %>% 
#   group_by(participant, block_id) %>%
#   mutate(
#     resp_hand = ifelse(key %in% c('f', 'g'), 'L', 'R'),
#     resp_fing = ifelse(key %in% c('f', 'j'), 'M', 'I'),
#     resp0_hand = lag(resp_hand),
#     resp0_fing = lag(resp_fing), 
#     hand_rep = resp_hand == resp0_hand,
#     fing_rep = resp_fing == resp0_fing,
#     both_rep = hand_rep & fing_rep,
#     H_rep_err = hand_rep & !both_rep & !accuracy,
#     F_rep_err = fing_rep & !both_rep & !accuracy,
#     HF_rep_err = both_rep & !accuracy, 
#     N_rep_err = !hand_rep & !fing_rep & !accuracy)
# plt_dat1 <- tmp_dat %>% mutate(X = PT, Y = accuracy) %>% smooth_obs(covar = c(covar_name), width = 100) %>% mutate(resp = 'correct')
# plt_dat2 <- tmp_dat %>% mutate(X = PT, Y = H_rep_err) %>% smooth_obs(covar = c(covar_name), width = 100) %>% mutate(resp = 'err, rep hand')
# plt_dat3 <- tmp_dat %>% mutate(X = PT, Y = F_rep_err) %>% smooth_obs(covar = c(covar_name), width = 100) %>% mutate(resp = 'err, rep finger')
# plt_dat4 <- tmp_dat %>% mutate(X = PT, Y = HF_rep_err) %>% smooth_obs(covar = c(covar_name), width = 100) %>% mutate(resp = 'err, rep both')
# plt_dat5 <- tmp_dat %>% mutate(X = PT, Y = N_rep_err) %>% smooth_obs(covar = c(covar_name), width = 100) %>% mutate(resp = 'err, rep none')

# plt_dat <- bind_rows(plt_dat1, plt_dat2, plt_dat3, plt_dat4, plt_dat5) %>% mutate(err0 = covar1)
  
# p <- plt_dat %>% 
#   ggplot(aes(x = PT, y = obs_mu, 
#     ymax = obs_mu + obs_se * 1.96, ymin = obs_mu - obs_se * 1.96, 
#     colour = resp, fill = resp, group = resp)) +
#   sw_geom + my_theme +  
#   scale_x_continuous("PT (ms)", breaks = seq(0, 1000, 100), labels = as.character(seq(0, 2000, 200))) +
#   scale_colour_manual(values = WONG[c(1, 2, 3, 5, 6)]) +
#   scale_fill_manual(values = WONG[c(1, 2, 3, 5, 6)]) +
#   facet_wrap(err0 ~ ., labeller = labeller(err0 = c('-0.5' = 'After correct', '0.5' = 'After error'))) +
#   theme(legend.position = 'right') +
#   ggsave(sprintf('fig/%s_resp_prop.jpg', exp_name), units = 'in', height = 2, width = 7.5, dpi = 1000)


# dat %>% 
#   mutate(err0 = factor(err0)) %>%
#   ggplot(aes(x = PT, fill = err0, colour = err0, group = err0)) +
#   geom_histogram(position = "identity", bins = 10, alpha = .5) +
#   scale_x_continuous("PT (ms)", breaks = seq(0, 1000, 100), labels = as.character(seq(0, 2000, 200))) +
#   coord_cartesian(xlim = c(0, 1000)) +
#   scale_colour_manual(values = WONG[c(4, 7)]) +
#   scale_fill_manual(values = WONG[c(4, 7)]) +
#   my_theme +
#  # hist_theme +
#   scale_y_continuous() +
#   ggsave(sprintf('fig/%s_PT_hist_err0.jpg', exp_name), units = 'in', height = 2, width = 3, dpi = 1000)


# dat %>%
#   filter(PT > 100) %>%
#   ggplot(aes(x = key_stim_rep, y = accuracy, colour = factor(err0), group = err0)) +
#   stat_summary(position = position_dodge(.1)) +
#   scale_colour_manual(values = WONG[c(4, 7)])

# dat %>%
#   filter(PT > 400) %>%
#   # group_by(key_stim_rep) %>%
#   # mutate(N = n()) %>%
#   # group_by(key_stim_rep, err0) %>%
#   # summarise(count = n(), prop = count / N) %>%
#   ggplot(aes(x = key_stim_rep, y = as.numeric(key_rep), colour = factor(err0), fill = factor(err0), group = err0)) +
#   stat_summary() +
#   #scale_y_continous("Probability of repeat") +
#   scale_x_discrete(labels = c('Should not repeat', 'Should Repeat')) +
#   scale_colour_manual(values = WONG[c(4, 7)], labels = c('After Correct', 'After Error')) +
#   scale_fill_manual(values = WONG[c(4, 7)], labels = c('After Correct', 'After Error'))

# # exp_name <- 'memory_err0'
# # fit1 <- read_rds(sprintf('fit/%s_rpm1_mu_1_sigma_1_beta_1_alpha_1.rds',exp_name))
# # fit_mix <- read_rds(sprintf('fit/%s_rpm_mix_mu_2_sigma_2_beta_1_alpha_1.rds',exp_name))

# # draws1 <- as_draws_df(fit1)
# # draws_mix <- as_draws_df(fit_mix)

# # plot_prep1(draws1, exp_name)

# # plot_prep_mix(draws_mix, exp_name)

# # get_lambda(draws_mix) %>% spread_draws(lambda) %>% median_qi()
# # get_mu(draws_mix, R_p = 2, R_mu = 2, covar = FALSE) %>% spread_draws(mu[r]) %>% median_qi()
# # get_sigma(draws_mix, R_p = 2, R_sigma = 2, covar = FALSE) %>% spread_draws(sigma[r]) %>% median_qi()



# # p <- mcmc_hist(fit_mix, pars = c('delta_lambda[1]'), binwidth = .005, 
# #   transform = list('delta_lambda[1]' = 'inv_logit')) + 
# #   my_theme + 
# #   geom_histogram(aes(x = inv_logit(rnorm(8e3, 2, .5))), alpha = .5, binwidth = .005) +
# #   scale_x_continuous(expression(lambda), breaks = seq(0, 1, .3)) +
# #   coord_cartesian(xlim = c(0, 1), ylim = c(0, 4000)) + 
# #   hist_theme +
# #   ggsave(sprintf('fig/%s_lambda_mcmc_hist_new_prior.jpg', mod_name), units = 'in', height = 1, width = 1.5, dpi = 1000)

