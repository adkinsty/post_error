library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(tidybayes) # Kay M (2020). tidybayes: Tidy Data and Geoms for Bayesian Models. doi: 10.5281/zenodo.1308151, R package version 2.3.1, http://mjskay.github.io/tidybayes/.
library(bayestestR)
library(bayesplot)
library(brms)
library(bayestestR)

covar_name <- 'err0'
mod_type <- 'rpm1_hier'

setwd('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/')
source('../../forced_response/func.R')

exp_name1 <- 'memory_color_short'
mod_name1 <- sprintf('%s_mu_%s_sigma_%s_beta_%s', mod_type, covar_name, covar_name, covar_name)
fit1 <- read_rds(sprintf('%s/fit/%s_%s.rds',exp_name1, exp_name1, mod_name1)) 
draws1 <- as_draws_df(fit1)
betas1 <- get_beta(draws1, R_beta = 1, covar = TRUE)
beta_diff1 <- get_diff(betas1, par = 'beta', R = 1) %>% subset_draws(variable = 'beta_diff')

exp_name2 <- 'memory_color_long'
mod_name2 <- sprintf('%s_mu_%s_sigma_%s_beta_%s', mod_type, covar_name, covar_name, covar_name)
fit2 <- read_rds(sprintf('%s/fit/%s_%s.rds',exp_name2, exp_name2, mod_name2)) 
draws2 <- as_draws_df(fit2)
betas2 <- get_beta(draws2, R_beta = 1, covar = TRUE)
beta_diff2 <- get_diff(betas2, par = 'beta', R = 1) %>% subset_draws(variable = 'beta_diff')

plt_dat <- tibble(
    diff = c(beta_diff1$beta_diff, beta_diff2$beta_diff),
    exp = rep(c('short', 'long'), each = nrow(beta_diff1)))

plt_dat %>%
    ggplot(aes(x = diff, group = exp, fill = exp, colour = exp)) +
    geom_histogram(
        binwidth = .005, 
        position = 'identity',
        alpha = .5) + 
    scale_x_continuous(expression(Delta[beta])) +
    scale_y_continuous('MCMC Draws') + 
    scale_fill_manual('ITI', values = c('red', 'blue3')) + 
    scale_colour_manual('ITI', values = c('red', 'blue3')) + 
    theme_classic(base_size = 20)

