# Constants and functions for analyses and plots of data from mot_simon studies
# Ordered alphabetically

library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(rstan)
library(loo) # Vehtari et al. (2020). “loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models.” R package version 2.4.1, https://mc-stan.org/loo/.
library(brms) # Bürkner P (2017). “brms: An R Package for Bayesian Multilevel Models Using Stan.” Journal of Statistical Software, 80(1), 1–28. doi: 10.18637/jss.v080.i01.
library(tidybayes) # Kay M (2020). tidybayes: Tidy Data and Geoms for Bayesian Models. doi: 10.5281/zenodo.1308151, R package version 2.3.1, http://mjskay.github.io/tidybayes/.
library(posterior) # Bürkner et al. (2020). “posterior: Tools for Working with Posterior Distributions.” R package version XXX, <URL: https://mc-stan.org/posterior>.

ci95 <- function(x) {
  return(quantile(x, prob = .95))
}


ci05 <- function(x) {
  return(quantile(x, prob = .05))
}


ci90 <- function(x) {
  return(quantile(x, prob = .9))
}


ci10 <- function(x) {
  return(quantile(x, prob = .1))
}


delta_summary <- function(fit, exp_name) {
  delta_mu <- fit %>% spread_draws(delta_mu[k, r]) %>% median_qi() %>% 
    write_csv(sprintf('tab/%s_delta_mu.csv',exp_name))
  delta_beta <- fit %>% spread_draws(delta_beta[k, r]) %>% median_qi() %>% 
    write_csv(sprintf('tab/%s_delta_beta.csv',exp_name))
  delta_sigma <- fit %>% spread_draws(delta_sigma[k, r]) %>% median_qi() %>% 
    write_csv(sprintf('tab/%s_delta_sigma.csv',exp_name))
  param <- list('mu' = delta_mu, 'sigma' = delta_sigma, 'beta' = delta_beta)
  return(param)
}


design_matrix <- function(x_names, dat) {
  # given dataset for one experiment and names of covariate columns
  # returns design matrix of predictors for input to rpm.stan 

  N <- nrow(dat)
  x <- matrix(rep(1, N), ncol = 1)

  for (name in x_names) {

    if (str_detect(name, '\\*')) {

      name <- str_replace_all(name, ' ', '')
      names <- str_split(name, '\\*')[[1]]
      new_x <- rep(1, N)

      for (tmp_name in names) {
        new_x <- new_x * pull(dat, tmp_name)
      }

    } else {

      new_x <- pull(dat, name)

    }

    x <- cbind(x, new_x)

  }

  return(x)

}


format_input <- function(dat,
  R_mu, R_sigma, R_beta, cov_mu, cov_sigma, cov_beta, cov_alpha, gq, prior_only) {

  # given dataset for one experiment and various model hyperparameters
  #   R_par controls number of processes for par
  #   cov_par controls whether to allow effects of covariates on par
  # return complete list of data formatted for input to rpm.stan

  input <- list()

  input$N <- nrow(dat)
  input$y <- dat$accuracy
  input$time <- dat$time
  input$err0 <- dat$err0

  input$R_mu <- R_mu
  input$R_sigma <- R_sigma
  input$R_beta <- R_beta

  input$x_mu <- design_matrix(cov_mu, dat)
  input$x_sigma <- design_matrix(cov_sigma, dat)
  input$x_beta <- design_matrix(cov_beta, dat)
  input$x_alpha <- design_matrix(cov_alpha, dat)

  input$K_mu <- ncol(input$x_mu)
  input$K_sigma <- ncol(input$x_sigma)
  input$K_beta <- ncol(input$x_beta)
  input$K_alpha <- ncol(input$x_alpha)

  input$gq <- gq
  input$prior_only <- prior_only

  return(input)
}

format_input_mix <- function(dat,
  R_mu, R_sigma, R_beta, cov_mu, cov_sigma, cov_beta, cov_alpha, cov_lambda, gq, prior_only) {

  # given dataset for one experiment and various model hyperparameters
  #   R_par controls number of processes for par
  #   cov_par controls whether to allow effects of covariates on par
  # return complete list of data formatted for input to rpm.stan

  input <- list()

  input$N <- nrow(dat)
  input$y <- dat$accuracy
  input$time <- dat$time

  input$R_mu <- R_mu
  input$R_sigma <- R_sigma
  input$R_beta <- R_beta

  input$x_mu <- design_matrix(cov_mu, dat)
  input$x_sigma <- design_matrix(cov_sigma, dat)
  input$x_beta <- design_matrix(cov_beta, dat)
  input$x_alpha <- design_matrix(cov_alpha, dat)
  input$x_lambda <- design_matrix(cov_lambda, dat)

  input$K_mu <- ncol(input$x_mu)
  input$K_sigma <- ncol(input$x_sigma)
  input$K_beta <- ncol(input$x_beta)
  input$K_alpha <- ncol(input$x_alpha)
  input$K_lambda <- ncol(input$x_lambda)

  input$gq <- gq
  input$prior_only <- prior_only

  return(input)
}

get_beta <- function(draws, R_beta = 1, covar = FALSE) {
  beta_draws <- draws %>% subset_draws(variable = 'delta_beta*', regex = TRUE)
  if (!covar) {
    beta_draws <- beta_draws %>%
      mutate_variables(`beta[h]` = inv_logit_scaled(`delta_beta[1,1]`))
    if (R_beta == 2) {
      beta_draws <- beta_draws %>%
        mutate_variables(`beta[g]` = inv_logit_scaled(`delta_beta[1,2]`))
    }
  } else {
    beta_draws <- beta_draws %>%
      mutate_variables(`beta[h,lo]` = inv_logit_scaled(`delta_beta[1,1]` + `delta_beta[2,1]` * -.5),
        `beta[h,hi]` = inv_logit_scaled(`delta_beta[1,1]` + `delta_beta[2,1]` * .5))
    if (R_beta == 2) {     
      beta_draws <- beta_draws %>% 
        mutate_variables(`beta[g,lo]` = inv_logit_scaled(`delta_beta[1,2]` + `delta_beta[2,2]` * -.5),
          `beta[g,hi]` = inv_logit_scaled(`delta_beta[1,2]` + `delta_beta[2,2]` * .5))
    }
  }
  return(beta_draws)
}

get_covariates <- function(exp_name) {
  covariates <- c('err0')
  return(covariates)
}

get_lambda <- function(draws) {
  lambda_draws <- draws %>% 
    subset_draws(variable = 'delta_lambda*', regex = TRUE) %>%
    mutate(lambda = inv_logit_scaled(`delta_lambda[1]`))
  return(lambda_draws)
}

get_mu <- function(draws, R_mu = 1, covar = TRUE) {
  mu_draws <- draws %>% subset_draws(variable = 'delta_mu*', regex = TRUE)
  if (!covar) {
    if (R_mu == 1) {
      mu_draws <- mu_draws %>%
      mutate_variables(mu = inv_logit_scaled(`delta_mu[1]`))
    } else {
      if (R_mu == 2) {
        mu_draws <- mu_draws %>%
          mutate_variables(`mu[h]` = inv_logit_scaled(`delta_mu[1,1]`),
            `mu[g]` = inv_logit_scaled(`delta_mu[1,2]`))
      }
    }
  } else {

    if (R_mu == 1) {
      mu_draws <- mu_draws %>%
        mutate_variables(`mu[hi]` = inv_logit_scaled(`delta_mu[1]` + `delta_mu[2]` * .5),
            `mu[lo]` = inv_logit_scaled(`delta_mu[1]` + `delta_mu[2]` * -.5))
    } else {
      if (R_mu == 2) {     
        mu_draws <- mu_draws %>% 
          mutate_variables(`mu[h,lo]` = inv_logit_scaled(`delta_mu[1,1]` + `delta_mu[2,1]` * -.5),
            `mu[h,hi]` = inv_logit_scaled(`delta_mu[1,1]` + `delta_mu[2,1]` * .5),
            `mu[g,lo]` = inv_logit_scaled(`delta_mu[1,2]` + `delta_mu[2,2]` * -.5),
            `mu[g,hi]` = inv_logit_scaled(`delta_mu[1,2]` + `delta_mu[2,2]` * .5))
      }
    }
  }
  return(mu_draws)
}

get_mu1 <- function(draws, covar = FALSE) {
  mu_draws <- draws %>% subset_draws(variable = 'delta_mu*', regex = TRUE)
  if (!covar) {
    mu_draws <- mu_draws %>%
      mutate_variables(`mu` = inv_logit_scaled(`delta_mu[1]`))
  } else {
    mu_draws <- mu_draws %>%
      mutate_variables(
        `mu[lo]` = inv_logit_scaled(`delta_mu[1]` + `delta_mu[2]` * -.5),
        `mu[hi]` = inv_logit_scaled(`delta_mu[1]` + `delta_mu[2]` * .5))
  }
  return(mu_draws)
}

get_mu_diff <- function(mu_draws, covar = TRUE) {
  if (!covar) {
    diff_draws <- mu_draws %>% mutate_variables(`mu_diff[lo]` = `mu[g]` - `mu[h]`)
  } else {
    diff_draws <- mu_draws %>%
      mutate_variables(
        `mu_diff[lo]` = `mu[g,lo]` - `mu[h,lo]`,
        `mu_diff[hi]` = `mu[g,hi]` - `mu[h,hi]`)
  }
  return(diff_draws)
}

get_mu_diff2 <- function(diff_draws) {
  diff_draws2 <- diff_draws %>% mutate_variables(`mu_diff2` = `mu_diff[hi]` - `mu_diff[lo]`)
  return(diff_draws2)
}

get_sigma <- function(draws, R_sigma = 1, covar = TRUE) {
  sigma_draws <- draws %>% subset_draws(variable = 'delta_sigma*', regex = TRUE)
  if (!covar) {
    if (R_sigma == 1) {
      sigma_draws <- sigma_draws %>%
      mutate_variables(sigma = inv_logit_scaled(`delta_sigma[1]`))
    } else {
      if (R_sigma == 2) {
        sigma_draws <- sigma_draws %>%
          mutate_variables(`sigma[h]` = inv_logit_scaled(`delta_sigma[1,1]`),
            `sigma[g]` = inv_logit_scaled(`delta_sigma[1,2]`))
      }
    }
  } else {

    if (R_sigma == 1) {
      sigma_draws <- sigma_draws %>%
        mutate_variables(`sigma[hi]` = inv_logit_scaled(`delta_sigma[1]` + `delta_sigma[2]` * .5),
            `sigma[lo]` = inv_logit_scaled(`delta_sigma[1]` + `delta_sigma[2]` * -.5))
    } else {
      if (R_sigma == 2) {     
        sigma_draws <- sigma_draws %>% 
          mutate_variables(`sigma[h,lo]` = inv_logit_scaled(`delta_sigma[1,1]` + `delta_sigma[2,1]` * -.5),
            `sigma[h,hi]` = inv_logit_scaled(`delta_sigma[1,1]` + `delta_sigma[2,1]` * .5),
            `sigma[g,lo]` = inv_logit_scaled(`delta_sigma[1,2]` + `delta_sigma[2,2]` * -.5),
            `sigma[g,hi]` = inv_logit_scaled(`delta_sigma[1,2]` + `delta_sigma[2,2]` * .5))
      }
    }
  }
  return(sigma_draws)
}

get_sigma1 <- function(draws, covar = FALSE) {
  sigma_draws <- draws %>% subset_draws(variable = 'delta_sigma*', regex = TRUE)
  if (!covar) {
    sigma_draws <- sigma_draws %>%
      mutate_variables(`sigma` = inv_logit_scaled(`delta_sigma[1]`))
  } else {
    sigma_draws <- sigma_draws %>%
      mutate_variables(
        `sigma[lo]` = inv_logit_scaled(`delta_sigma[1]` + `delta_sigma[2]` * -.5),
        `sigma[hi]` = inv_logit_scaled(`delta_sigma[1]` + `delta_sigma[2]` * .5))
  }
  return(sigma_draws)
}


inv_logit <- function(x) {
  y <- exp(x) / (1 + exp(x))
  return(y)
}

multiverse <- function(dat = c(), cov_param = c('mu', 'sigma', 'beta'), R_param = c('sigma', 'beta')) {
  # given dataset for many experiments
  # re-run model under many hyperparameter settings
  # concatenate rpm outputs
  # return tibble of hyperparameters and loo model outputs

  cov_param <- rje::powerSet(cov_param) # possible param to apply covariates
  R_param <- rje::powerSet(R_param) # possible param to use more than 1 proccess
  all_exp <- unique(dat$exp)

  model_fits <- tibble()

  for (exp_name in all_exp) {

    setting <- 1
    rpm_arg = list()

    exp_dat <- dat %>% filter(exp == exp_name)
    covariates <- get_covariates(exp_name)

    for (i in 1:length(cov_param)) { 
      for (j in 1:length(R_param)) {
        print('**********************************************************************************************')
        print('**********************************************************************************************')
        print(paste0('Experiment: ',exp_name))
        print(paste0('Permutation:  ', setting))

        rpm_arg[['dat']] = exp_dat
        rpm_arg[['cov_mu']] = covariates
        for (par_name in cov_param[[i]]) {
          rpm_arg[[paste0('cov_',par_name)]] = covariates}
        for (par_name in R_param[[j]]) {
          rpm_arg[[paste0('R_',par_name)]] = 2}

        rpm_output <- rlang::invoke(rpm, rpm_arg)

        model_fits <- model_fits %>% rbind(rpm_output)


        setting <- setting + 1
        rpm_arg <- list()
        rpm_output <- c()
        gc()
      }
    }
  }
  write_rds(model_fits, 'loo/model_fits.rds')
  return(model_fits)
}

plot_beta <- function(draws, exp_name) {

  betas <- get_beta(draws)

  p <- betas %>% 
    spread_draws(beta[r]) %>%
    ggplot(aes(x = beta, fill = r, colour = r)) + 
    stat_halfeye(point_size = 1, slab_alpha = .5) +
    scale_fill_manual('Response', values = WONG[c(4, 8)]) +
    scale_colour_manual('Response', values = WONG[c(4, 8)]) +
    scale_x_continuous(expression(paste('Expression Probability (',beta,')')), breaks = seq(0, 1, .05)) +  
    coord_cartesian(xlim = c(.7, 1)) + 
    my_theme +  
    theme(axis.line.y = element_blank(), 
      axis.text.y = element_blank(), 
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank())
      
  path <- sprintf('fig/%s_beta.jpg', exp_name)
  ggsave(plot = p, filename = path, units = 'in', height = 1.25, width = 2.75)
  
  return(p)
}

plot_mu <- function(draws, exp_name) {

  mus <- get_mu(draws)

  p <- mus %>% 
    spread_draws(mu[r, v]) %>% mutate(v = factor(v, c('lo', 'hi'))) %>%
    median_qi() %>% ungroup() %>%
    mutate(y = 1:4) %>%
    ggplot(aes(x = mu*1000, xmax = .upper*1000, xmin = .lower*1000, y = y, alpha = v, fill = r, colour = r)) + 
    geom_pointinterval(size = .5) +
    scale_fill_manual('Response', values = WONG[c(4, 8)]) +
    scale_colour_manual('Response', values = WONG[c(4, 8)]) +
    scale_alpha_manual('', values = c(.5, 1)) +
    scale_y_continuous("Parameter", expand = c(.2, .2), breaks = 1:4, labels = c(expression(mu["g,no"]), expression(mu["g,yes"]), expression(mu["h,no"]), expression(mu["h,yes"]))) +
    scale_x_continuous("Time (ms)", breaks = seq(0, 1000, 100)) +  
    coord_cartesian(xlim = c(0, 1000)) + 
    my_theme +  
    theme(panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank())
      
  path <- sprintf('fig/%s_mu.jpg',exp_name)
  ggsave(plot = p, filename = path, units = 'in', height = 1.5, width = 3, dpi = 1000)

  return(p)
}

plot_prep1 <- function(draws, exp_name) {

  mus <- draws %>% subset_draws(variable = 'delta_mu*', regex = TRUE) %>% mutate(mu = inv_logit_scaled(`delta_mu[1]`))
  sigmas <- draws %>% subset_draws(variable = 'delta_sigma*', regex = TRUE) %>% mutate(sigma = inv_logit_scaled(`delta_sigma[1]`))

  prep <- mus %>% 
    inner_join(sigmas) %>% 
    median_qi() %>%
    summarise(t = seq(0, 1, .001),
      dens = dnorm(t, mu, sigma))

  p <- prep %>% 
    ggplot(aes(x = t*1000, y = dens)) + 
    geom_line() +
    scale_x_continuous('Time (ms)', breaks = seq(0, 1000, 100)) +  
    coord_cartesian(xlim = c(0, 1000), ylim = c(0, 5)) + 
    my_theme +  
    theme(axis.line.y = element_blank(), 
      panel.grid.minor.y = element_blank(),
      panel.grid.major = element_blank())

  path <- sprintf('fig/%s_prep_dens1.jpg', exp_name)
  ggsave(plot = p, filename = path, units = 'in', height = 1.5, width = 3, dpi = 1000)

  return(p)
}

plot_prep1_cov <- function(draws, exp_name) {

  mus <- draws %>% subset_draws(variable = 'delta_mu*', regex = TRUE) %>% 
    mutate(`mu[hi]` = inv_logit_scaled(`delta_mu[1]` + `delta_mu[2]` * .5),
      `mu[lo]` = inv_logit_scaled(`delta_mu[1]` + `delta_mu[2]` * -.5)) %>% spread_draws(mu[v])

  sigmas <- draws %>% subset_draws(variable = 'delta_sigma*', regex = TRUE) %>% 
    mutate(`sigma[hi]` = inv_logit_scaled(`delta_sigma[1]` + `delta_sigma[2]` * .5),
      `sigma[lo]` = inv_logit_scaled(`delta_sigma[1]` + `delta_sigma[2]` * -.5)) %>% spread_draws(sigma[v])

  prep <- mus %>% 
    inner_join(sigmas) %>% 
    median_qi() %>%
    group_by(v) %>%
    summarise(t = seq(0, 1, .001),
      dens = dnorm(t, mu, sigma))

  p <- prep %>% 
    ggplot(aes(x = t*1000, y = dens, fill = v, colour = v, group = v)) + 
    geom_line() +
    scale_fill_manual('Err0', values = WONG[c(4, 7)]) +
    scale_colour_manual('Err0', values = WONG[c(4, 7)]) +
    scale_x_continuous('Time (ms)', breaks = seq(0, 1000, 100)) +  
    coord_cartesian(xlim = c(0, 1000), ylim = c(0, 5)) + 
    my_theme +  
    theme(axis.line.y = element_blank(), 
      panel.grid.minor.y = element_blank(),
      panel.grid.major = element_blank())

  path <- sprintf('fig/%s_prep_dens1_cov.jpg', exp_name)
  ggsave(plot = p, filename = path, units = 'in', height = 1.5, width = 3, dpi = 1000)

  return(p)
}

plot_prep2_cov <- function(draws, exp_name) {
  sigmas <- get_sigma(draws) %>% spread_draws(sigma[r, v])
  mus <- get_mu(draws) %>% spread_draws(mu[r, v])
  prep <- mus %>% 
    inner_join(sigmas) %>% 
    median_qi() %>%
    group_by(r, v) %>%
    summarise(t = seq(0, 1, .001),
      dens = dnorm(t, mu, sigma))

  p <- prep %>% mutate(v = factor(v, c('lo', 'hi'))) %>%
    ggplot(aes(x = t*1000, y = dens, fill = v, colour = v, group = v)) + 
    geom_line() +
    scale_fill_manual('Response', values = WONG[c(4, 7)]) +
    scale_colour_manual('Response', values = WONG[c(4, 7)]) +
    scale_alpha_manual('Reward', values = c(.5, 1)) +
    scale_x_continuous('Time (ms)', breaks = seq(0, 1000, 100)) +  
    coord_cartesian(xlim = c(0, 1000), ylim = c(0, 10)) + 
    my_theme +  
    theme(axis.line.y = element_blank(), 
      panel.grid.minor.y = element_blank(),
      panel.grid.major = element_blank())

  path <- sprintf('fig/%s_prep_dens2_cov.jpg', exp_name)
  ggsave(plot = p, filename = path, units = 'in', height = 1.5, width = 3, dpi = 1000)

  return(p)
}

plot_prep_mix <- function(draws, exp_name) {
  mus <- draws %>% subset_draws(variable = 'delta_mu*', regex = TRUE) %>% 
    mutate(mu1 = inv_logit_scaled(`delta_mu[1,1]`),
      mu2 = inv_logit_scaled(`delta_mu[1,2]`))
  sigmas <- draws %>% subset_draws(variable = 'delta_sigma*', regex = TRUE) %>% 
    mutate(sigma1 = inv_logit_scaled(`delta_sigma[1,1]`),
      sigma2 = inv_logit_scaled(`delta_sigma[1,2]`))
  lambdas <- draws %>% subset_draws(variable = 'delta_lambda*', regex = TRUE) %>% 
    mutate(lambda = inv_logit_scaled(`delta_lambda[1]`))
  
  prep <- mus %>% 
    inner_join(sigmas) %>% 
    inner_join(lambdas) %>%
    median_qi() %>%
    summarise(t = seq(0, 1, .001),
      dens = lambda * dnorm(t, mu1, sigma1) + (1 - lambda) * dnorm(t, mu2, sigma2), 
      dens1 = lambda * dnorm(t, mu1, sigma1),
      dens2 = (1 - lambda) * dnorm(t, mu2, sigma2))

  p <- prep %>%
    ggplot(aes(x = t*1000)) + 
    geom_line(aes(y = dens), colour = 'green') +
    geom_line(aes(y = dens1), colour = 'blue') +
    geom_line(aes(y = dens2), colour = 'red') +
    scale_x_continuous('Time (ms)', breaks = seq(0, 1000, 100)) +  
    coord_cartesian(xlim = c(0, 1000), ylim = c(0, 5)) + 
    my_theme +  
    theme(axis.line.y = element_blank(), 
      panel.grid.minor.y = element_blank(),
      panel.grid.major = element_blank())

  path <- sprintf('fig/%s_prep_dens_mix.jpg', exp_name)
  ggsave(plot = p, filename = path, units = 'in', height = 1.5, width = 3, dpi = 1000)

  return(p)
}

plot_sigma <- function(draws, exp_name) {

  sigmas <- get_sigma(draws)

  p <- sigmas %>% 
    spread_draws(sigma[r, v]) %>% mutate(v = factor(v, c('lo', 'hi'))) %>%
    median_qi() %>% ungroup() %>%
    mutate(y = 1:4) %>%
    ggplot(aes(x = sigma*1000, xmax = .upper*1000, xmin = .lower*1000, y = y, alpha = v, fill = r, colour = r)) + 
    geom_pointinterval(size = .5) +
    scale_fill_manual('Response', values = WONG[c(4, 8)]) +
    scale_colour_manual('Response', values = WONG[c(4, 8)]) +
    scale_alpha_manual('', values = c(.5, 1)) +
    scale_y_continuous("Parameter", expand = c(.2, .2), breaks = 1:4, labels = c(expression(sigma["g,no"]), expression(sigma["g,yes"]), expression(sigma["h,no"]), expression(sigma["h,yes"]))) +
    scale_x_continuous("Time (ms)", breaks = seq(0, 700, 100)) +  
    coord_cartesian(xlim = c(0, 700)) + 
    my_theme +  
    theme(panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank())
      
  path <- sprintf('fig/%s_sigma.jpg',exp_name)
  ggsave(plot = p, filename = path, units = 'in', height = 1.5, width = 2, dpi = 1000)

  return(p)
}

reflect <- function(x, goals) {
  # convert responses to accuracies
  x = ifelse(goals == 0, 1 - x, x)
  return(x)
}

rpm1 <- function(dat,
  exp_name = NA,
  R_mu = 1, R_sigma = 1, R_beta = 1, 
  cov_mu = c(), cov_sigma = c(), cov_beta = c(), cov_alpha = c(), 
  gq = 1, prior_only = 0, re_run = FALSE,
  chains = 4, cores = 2, iter = 3000, warmup = 1000) {

  # R interface for the single process response preparation model in rpm1.stan
  # Given dataset and settings for hyperparameters
  # Do Sample from model and do LOO evaluation and write results
  # Return data frame row with hyperparam. and loo fit output

  if (is.na(exp_name)) {
    exp_name <- dat$exp[1]
  }  
  
  fit_filename <- sprintf('fit/rpm1_%s_mu_1%s_sigma_1%s_beta_1%s_alpha_1%s.rds', 
    exp_name, str_C(cov_mu), str_C(cov_sigma), str_C(cov_beta), str_C(cov_alpha))
  
  loo_filename <- sprintf('loo/loo_rpm1_%s_mu_1%s_sigma_1%s_beta_1%s_alpha_1%s.rds', 
    exp_name, str_C(cov_mu), str_C(cov_sigma), str_C(cov_beta), str_C(cov_alpha))

  if (re_run | !file.exists(fit_filename)) {

    print(paste0(fit_filename, ' does not exist yet.'))
    print('Compiling...')
    
    stan_code <- stan_model('rpm1.stan')

    stan_dat <- format_input(
      dat, 
      R_mu, R_sigma, R_beta, 
      cov_mu, cov_sigma, cov_beta, cov_alpha, 
      gq, prior_only)

    print('Sampling...')
    stan_fit <- sampling(
      object = stan_code, 
      dat = stan_dat, 
      chains = chains, 
      cores = cores, 
      iter = iter, 
      warmup = warmup)

    print('Writing stanfit...')
    write_rds(stan_fit, fit_filename)

    print('Evaluating loo...')
    loo_fit <- loo(stan_fit, cores = 4)

    print('Writing loo...')
    write_rds(loo_fit, loo_filename)

  
  } else {
    if (re_run | !file.exists(loo_filename)) {

      print(paste0(loo_filename, ' does not exist yet.'))
      print('Reading stanfit...')
      stan_fit <- read_rds(fit_filename)

      print('Evaluating loo...')
      loo_fit <- loo(stan_fit, cores = 4)

      print('Writing loo...')
      write_rds(loo_fit, loo_filename)

    } else {

      print(paste0(loo_filename, ' already exists!'))
      print('Reading loo...')
      loo_fit <- read_rds(loo_filename)

    }
  }  
  
  output <- data.frame() %>%
    do(loo = loo_fit) %>%
    mutate(
      exp = exp_name,
      R_mu, R_sigma, R_beta,
      cov_mu = str_C(cov_mu), cov_sigma = str_C(cov_sigma),  
      cov_beta = str_C(cov_beta), cov_alpha = str_C(cov_alpha),
      filename = fit_filename)

  return(output)
}

rpm2 <- function(dat,
  exp_name = NA,
  R_mu = 2, R_sigma = 2, R_beta = 2, 
  cov_mu = c(), cov_sigma = c(), cov_beta = c(), cov_alpha = c(), 
  gq = 1, prior_only = 0, re_run = FALSE,
  chains = 4, cores = 2, iter = 3000, warmup = 1000) {

  # R interface for the dual process response preparation model in rpm2.stan
  # Given dataset and settings for hyperparameters
  # Do Sample from model and do LOO evaluation and write results
  # Return data frame row with hyperparam. and loo fit output

  if (is.na(exp_name)) {
    exp_name <- dat$exp[1]
  }  
  fit_filename <- sprintf('fit/rpm2_%s_mu_%s%s_sigma_%s%s_beta_%s%s_alpha_1%s.rds', 
    exp_name, R_mu, str_C(cov_mu), R_sigma, str_C(cov_sigma), R_beta, str_C(cov_beta), str_C(cov_alpha))
  
  loo_filename <- sprintf('loo/loo_rpm2_%s_mu_%s%s_sigma_%s%s_beta_%s%s_alpha_1%s.rds', 
    exp_name, R_mu, str_C(cov_mu), R_sigma, str_C(cov_sigma), R_beta, str_C(cov_beta), str_C(cov_alpha))

  if (re_run | !file.exists(fit_filename)) {

    print(paste0(fit_filename, ' does not exist yet.'))
    print('Compiling...')
    
    stan_code <- stan_model('rpm2.stan')

    stan_dat <- format_input(
      dat, 
      R_mu, R_sigma, R_beta, 
      cov_mu, cov_sigma, cov_beta, cov_alpha, 
      gq, prior_only)

    print('Sampling...')
    stan_fit <- sampling(
      object = stan_code, 
      dat = stan_dat, 
      chains = chains, 
      cores = cores, 
      iter = iter, 
      warmup = warmup)

    print('Writing stanfit...')
    write_rds(stan_fit, fit_filename)

    print('Evaluating loo...')
    loo_fit <- loo(stan_fit, cores = 4)

    print('Writing loo...')
    write_rds(loo_fit, loo_filename)

  
  } else {
    if (re_run | !file.exists(loo_filename)) {

      print(paste0(loo_filename, ' does not exist yet.'))
      print('Reading stanfit...')
      stan_fit <- read_rds(fit_filename)

      print('Evaluating loo...')
      loo_fit <- loo(stan_fit, cores = 4)

      print('Writing loo...')
      write_rds(loo_fit, loo_filename)

    } else {

      print(paste0(loo_filename, ' already exists!'))
      print('Reading loo...')
      loo_fit <- read_rds(loo_filename)

    }
  }  
  
  output <- data.frame() %>%
    do(loo = loo_fit) %>%
    mutate(
      exp = exp_name,
      R_mu, R_sigma, R_beta,
      cov_mu = str_C(cov_mu), cov_sigma = str_C(cov_sigma),  
      cov_beta = str_C(cov_beta), cov_alpha = str_C(cov_alpha),
      filename = fit_filename)

  return(output)
}

rpm_ctrl <- function(dat,
  exp_name = NA,
  R_mu = 1, R_sigma = 1, R_beta = 1, 
  cov_mu = c(), cov_sigma = c(), cov_beta = c(), cov_alpha = c(), 
  gq = 1, prior_only = 0, re_run = FALSE,
  chains = 4, cores = 2, iter = 3000, warmup = 1000) {

  # R interface for the response preparation model in rpm_ctrl.stan
  # Given dataset and settings for hyperparameters
  # Do Sample from model and do LOO evaluation and write results
  # Return data frame row with hyperparam. and loo fit output

  if (is.na(exp_name)) {
    exp_name <- dat$exp[1]
  }  
  fit_filename <- sprintf('fit/rpm_ctrl_%s_mu_%s%s_sigma_%s%s_beta_%s%s_alpha_1%s.rds', 
    exp_name, R_mu, str_C(cov_mu), R_sigma, str_C(cov_sigma), R_beta, str_C(cov_beta), str_C(cov_alpha))
  
  loo_filename <- sprintf('loo/loo_rpm_ctrl_%s_mu_%s%s_sigma_%s%s_beta_%s%s_alpha_1%s.rds', 
    exp_name, R_mu, str_C(cov_mu), R_sigma, str_C(cov_sigma), R_beta, str_C(cov_beta), str_C(cov_alpha))

  if (re_run | !file.exists(fit_filename)) {

    print(paste0(fit_filename, ' does not exist yet.'))
    print('Compiling...')
    
    stan_code <- stan_model('rpm_ctrl.stan')

    stan_dat <- format_input(
      dat, 
      R_mu, R_sigma, R_beta, 
      cov_mu, cov_sigma, cov_beta, cov_alpha, 
      gq, prior_only)

    print('Sampling...')
    stan_fit <- sampling(
      object = stan_code, 
      dat = stan_dat, 
      chains = chains, 
      cores = cores, 
      iter = iter, 
      warmup = warmup)

    print('Writing stanfit...')
    write_rds(stan_fit, fit_filename)

    print('Evaluating loo...')
    loo_fit <- loo(stan_fit, cores = 4)

    print('Writing loo...')
    write_rds(loo_fit, loo_filename)

  
  } else {
    if (re_run | !file.exists(loo_filename)) {

      print(paste0(loo_filename, ' does not exist yet.'))
      print('Reading stanfit...')
      stan_fit <- read_rds(fit_filename)

      print('Evaluating loo...')
      loo_fit <- loo(stan_fit, cores = 4)

      print('Writing loo...')
      write_rds(loo_fit, loo_filename)

    } else {

      print(paste0(loo_filename, ' already exists!'))
      print('Reading loo...')
      loo_fit <- read_rds(loo_filename)

    }
  }  
  
  output <- data.frame() %>%
    do(loo = loo_fit) %>%
    mutate(
      exp = exp_name,
      R_mu, R_sigma, R_beta,
      cov_mu = str_C(cov_mu), cov_sigma = str_C(cov_sigma),  
      cov_beta = str_C(cov_beta), cov_alpha = str_C(cov_alpha),
      filename = fit_filename)

  return(output)
}

rpm_mix <- function(dat,
  exp_name = NA,
  R_mu = 2, R_sigma = 2, R_beta = 1,
  cov_mu = c(), cov_sigma = c(), cov_beta = c(), cov_alpha = c(), cov_lambda = c(),
  gq = 1, prior_only = 0, re_run = FALSE,
  chains = 4, cores = 2, iter = 3000, warmup = 1000) {

  # R interface for the response preparation model in rpm.stan
  # Given dataset and settings for hyperparameters
  # Do Sample from model and do LOO evaluation and write results
  # Return data frame row with hyperparam. and loo fit output

  if (is.na(exp_name)) {
    exp_name <- dat$exp[1]
  }

  fit_filename <- sprintf('fit/rpm_mix_%s_mu_%s%s_sigma_%s%s_beta_1%s_alpha_1%s_lambda_1%s.rds', 
    exp_name, R_mu, str_C(cov_mu), R_sigma, str_C(cov_sigma), str_C(cov_beta), str_C(cov_alpha), str_C(cov_lambda))
  
  loo_filename <- sprintf('loo/loo_rpm_mix_%s_mu_%s%s_sigma_%s%s_beta_1%s_alpha_1%s_lambda_1%s.rds', 
    exp_name, R_mu, str_C(cov_mu), R_sigma, str_C(cov_sigma), str_C(cov_beta), str_C(cov_alpha), str_C(cov_lambda))

  if (re_run | !file.exists(fit_filename)) {

    print(paste0(fit_filename, ' does not exist yet.'))
    print('Compiling...')
    
    stan_code <- stan_model('rpm_mix.stan')

    stan_dat <- format_input(
      dat, 
      R_mu, R_sigma, R_beta, 
      cov_mu, cov_sigma, cov_beta, cov_alpha, cov_lambda,
      gq, prior_only)

    print('Sampling...')
    stan_fit <- sampling(
      object = stan_code, 
      dat = stan_dat, 
      chains = chains, 
      cores = cores, 
      iter = iter, 
      warmup = warmup)

    print('Writing stanfit...')
    write_rds(stan_fit, fit_filename)

    print('Evaluating loo...')
    loo_fit <- loo(stan_fit, cores = 4)

    print('Writing loo...')
    write_rds(loo_fit, loo_filename)

  
  } else {
    if (re_run | !file.exists(loo_filename)) {

      print(paste0(loo_filename, ' does not exist yet.'))
      print('Reading stanfit...')
      stan_fit <- read_rds(fit_filename)

      print('Evaluating loo...')
      loo_fit <- loo(stan_fit, cores = 4)

      print('Writing loo...')
      write_rds(loo_fit, loo_filename)

    } else {

      print(paste0(loo_filename, ' already exists!'))
      print('Reading loo...')
      loo_fit <- read_rds(loo_filename)

    }
  }  
  
  output <- data.frame() %>%
    do(loo = loo_fit) %>%
    mutate(
      exp = exp_name,
      R_mu, R_sigma, R_beta,
      cov_mu = str_C(cov_mu), cov_sigma = str_C(cov_sigma),  
      cov_beta = str_C(cov_beta), cov_alpha = str_C(cov_alpha),
      cov_lambda = str_C(cov_lambda),
      filename = fit_filename)

  return(output)
}

se <- function(x, na.rm = TRUE) {
  sd(x, na.rm = TRUE) / sqrt(length(x))
}

smooth_obs <- function(dat, lb = 0, ub = 1000, cov1 = NA, cov2 = NA) {
  # given dataset and model for a single experiment
  # run sliding window analysis on accuracy over time
  # return smoothed time-series of accuracy
  
  obs_pred <- dat

  exp_name <- dat$exp[1]

  dat$cov1 <- dat[cov1]

  if (is.na(cov2)) {

    tmp <- obs_pred %>% group_by(err0)
  } else {
    if (third_var == 'err00') {
      tmp <- obs_pred %>% group_by(err0, err00)
    } else {
      if (third_var == )
      print('group_by undefined for this third variable')
    }
  }

  sw_dat <- tmp %>% summarise(
    t = lb:ub, 
    obs_mu = sw_smooth(x = PT, y = accuracy),
    obs_se = sw_smooth(x = PT, y = accuracy, FUN = se))

  return(sw_dat)
}

smooth_obs_pred <- function(dat, fit, lb = 0, ub = 1000) {
  # given dataset and model for a single experiment
  # run sliding window analysis on accuracy over time
  # return smoothed time-series of accuracy
  
  pred <- t(rstan::extract(fit, pars = 'theta')$theta)

  obs_pred <- dat %>% 
    mutate(
      theta = apply(X = pred, FUN = mean, MARGIN = 1),
      high = apply(X = pred, FUN = ci95, MARGIN = 1),
      low = apply(X = pred, FUN = ci05, MARGIN = 1))


  exp_name <- dat$exp[1]

  tmp <- obs_pred %>% group_by(err0)

  sw_dat <- tmp %>% summarise(
    t = lb:ub, 
    obs_mu = sw_smooth(x = PT, y = accuracy),
    obs_se = sw_smooth(x = PT, y = accuracy, FUN = se),
    pred  = sw_smooth(x = PT, y = theta),
    upper = sw_smooth(x = PT, y = high),
    lower = sw_smooth(x = PT, y = low))

  return(sw_dat)
}


str_C <- function(v, sep = '+') {
  # given a vector of strings
  # return a single string with elements separated by sep
  return(str_c(v, collapse = sep))
}


sw_smooth <- function(x, y, lb = 0, ub = 1000, width = 100, FUN = mean) {
  # performing sliding window smoothing
  sw <- c()

  for (i in lb:ub) {

    lower <- i - (width / 2) 
    upper <- i + (width / 2)

    win <- y[which(x <= upper & x >= lower)]

    if (!is_empty(sample)) {
      sw <- append(sw, FUN(win, na.rm = TRUE))
    } else {
      sw <- append(sw, NaN) 
    }

  }
  return(sw)
}


to_name <- function(text, arg) {
  ctext <- sprintf(text, arg)
  name <-  as.name(ctext)
  return(name)
}

WONG <- c('#000000', '#e69f00', '#56b4e9', '#009e73', '#f0e442', '#0072b2', '#d55e00', '#cc79a7')

my_theme <- list(
  # re-usable theme for all plots
  theme_minimal(base_size = 10, base_family = 'sans'),
  theme(
    legend.position = 'none',
    legend.key.size = unit(.4, 'cm'), 
    axis.ticks.length = unit(0.1, 'cm'), 
    axis.ticks = element_line(size = .5),
    axis.line = element_line(size = .5),
    panel.grid.minor = element_blank()))


sw_geom <- list(
  # re-usable geoms for smoothed accuracy plots
    geom_line(size = .25),
    geom_ribbon(
      colour = NA, alpha = .1),
    scale_x_continuous(
      'PT (ms)', 
      breaks = seq(0, 1000, 100)),
    scale_y_continuous(
      'Accuracy', 
      breaks = seq(0, 1, .25), 
      labels = c('0', '.25', '.5', '.75', '1')),
    scale_colour_manual(
      values = c(WONG[4], WONG[7])),
    scale_fill_manual(
      values = c(WONG[4], WONG[7])),
    scale_alpha_manual(
      values = c(.5, 1)),
    coord_cartesian(
      xlim = c(0, 1000), 
      ylim = c(0, 1)))
