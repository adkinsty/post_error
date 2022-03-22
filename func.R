# Constants and functions for analyses and plots of data from forced response studies
# Ordered alphabetically

library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(rstan)
library(loo) # Vehtari et al. (2020). “loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models.” R package version 2.4.1, https://mc-stan.org/loo/.
library(brms) # Bürkner P (2017). “brms: An R Package for Bayesian Multilevel Models Using Stan.” Journal of Statistical Software, 80(1), 1–28. doi: 10.18637/jss.v080.i01.
library(tidybayes) # Kay M (2020). tidybayes: Tidy Data and Geoms for Bayesian Models. doi: 10.5281/zenodo.1308151, R package version 2.3.1, http://mjskay.github.io/tidybayes/.
library(posterior) # Bürkner et al. (2020). “posterior: Tools for Working with Posterior Distributions.” R package version XXX, <URL: https://mc-stan.org/posterior>.
library(bayesplot)

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

fit_rpm <- function(dat,
  exp_name = NA,
  mod_name = NA,
  covar_name = NA,
  gq = 1, prior_only = 0, re_run = FALSE,
  chains = 4, cores = 1, iter = 4000, warmup = 2000) {

  if (mod_name == 'rpm2_mu_cov') {
    fit_filename <- sprintf('fit/%s_%s_mu_%s_sigma_0_beta_0.rds', exp_name, mod_name, covar_name)
    loo_filename <- sprintf('loo/loo_%s_%s_mu_%s_sigma_0_beta_0.rds', exp_name, mod_name, covar_name)
  } else {
    fit_filename <- sprintf('fit/%s_%s_mu_%s_sigma_%s_beta_%s.rds', exp_name, mod_name, covar_name, covar_name, covar_name)
    loo_filename <- sprintf('loo/loo_%s_%s_mu_%s_sigma_%s_beta_%s.rds', exp_name, mod_name, covar_name, covar_name, covar_name)
  }


  if (re_run | !file.exists(fit_filename)) {

    print(paste0(fit_filename, ' does not exist yet.'))
    print('Compiling...')
    
    stan_code <- stan_model(sprintf('../../../forced_response/rpm/%s.stan', mod_name))

    stan_dat <- format_input(dat, covar_name, gq, prior_only)

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
    loo_fit <- loo(stan_fit, cores = 1)

    print('Writing loo...')
    write_rds(loo_fit, loo_filename)

  
  } else {
    if (re_run | !file.exists(loo_filename)) {

      print(paste0(loo_filename, ' does not exist yet.'))
      print('Reading stanfit...')
      stan_fit <- read_rds(fit_filename)

      print('Evaluating loo...')
      loo_fit <- loo(stan_fit, cores = 1)

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
      filename = fit_filename)

  return(output)
}

format_input <- function(dat, covar_name, gq, prior_only) {

  # temporary function for setting up input to hierarchical (varying intercepts) rpm w/ post-error effects

  input <- list()

  input$N <- nrow(dat)
  input$J <- length(unique(dat$participant))
  input$y <- dat$y
  input$time <- dat$time
  input$id <- as.numeric(as.factor(dat$participant))
  input$gq <- gq
  input$prior_only <- prior_only

  if (!str_detect(dat$exp[1], 'memory')) {
    # temporary hack to infer whether there is conflict/congruency
    input$congr <- dat$congr
  }
  if (dat$exp[1] %in% c('color')) {
    input$congr_rep <- dat$congr_rep
  }

  if (covar_name == 'err0') {
    input$x <- dat$err0
  }
  if (covar_name == 'rew') {
    input$x <- dat$rew
  }
  if (covar_name == 'congr0') {
    input$x <- dat$congr0
  }
  if (covar_name == 'congr') {
    input$x <- dat$congr
  }
  if (covar_name == 'err00') {
    input$x <- dat$err00
  }

  return(input)
}

get_beta <- function(draws, R_beta = 2, covar = TRUE, transform = inv_logit) {
  # variability in preparation demands
  raw_draws <- draws #%>% subset_draws(variable = 'delta_beta*', regex = TRUE)

  if (!covar) {
    if (R_beta == 1) {
      beta_draws <- raw_draws %>% mutate_variables(beta = transform(beta_loc))
    } else {
      beta_draws <- raw_draws %>%
        mutate_variables(
          `beta[h]` = transform(`beta_loc[1]`),
          `beta[g]` = transform(`beta_loc[2]`))
    }
  } else {
    if (R_beta == 1) {
      beta_draws <- raw_draws %>%
        mutate_variables(
          `beta[lo]` = transform(beta_loc + delta_beta_loc * -.5),
          `beta[hi]` = transform(beta_loc + delta_beta_loc * .5))
    } else {
      beta_draws <- raw_draws %>% 
        mutate_variables(
          `beta[h,lo]` = transform(`beta_loc[1]` + `delta_beta_loc[1]` * -.5),
          `beta[h,hi]` = transform(`beta_loc[1]` + `delta_beta_loc[1]` * .5),
          `beta[g,lo]` = transform(`beta_loc[2]` + `delta_beta_loc[2]` * -.5),
          `beta[g,hi]` = transform(`beta_loc[2]` + `delta_beta_loc[2]` * .5))
    }
  }
  return(beta_draws)
}

get_diff <- function(par_draws, par = 'mu', R = 2) {
  # get covariate effects on parameters
  if (par == 'mu') {
    if (R == 2) {
      diff_draws <- par_draws %>%
        mutate_variables(
          `mu_diff[h]` = `mu[h,hi]` - `mu[h,lo]`,
          `mu_diff[g]` = `mu[g,hi]` - `mu[g,lo]`)
    } else {
      diff_draws <- par_draws %>%
        mutate_variables(
          `mu_diff` = `mu[hi]` - `mu[lo]`)
    }
  }
  if (par == 'sigma') {
    if (R == 2) {
      diff_draws <- par_draws %>%
        mutate_variables(
          `sigma_diff[h]` = `sigma[h,hi]` - `sigma[h,lo]`,
          `sigma_diff[g]` = `sigma[g,hi]` - `sigma[g,lo]`)
    } else {
      diff_draws <- par_draws %>%
        mutate_variables(
          `sigma_diff` = `sigma[hi]` - `sigma[lo]`)
    }
  }
  if (par == 'beta') {
    if (R == 2) {
      diff_draws <- par_draws %>%
        mutate_variables(
          `beta_diff[h]` = `beta[h,hi]` - `beta[h,lo]`,
          `beta_diff[g]` = `beta[g,hi]` - `beta[g,lo]`)
    } else {
      diff_draws <- par_draws %>%
        mutate_variables(
          `beta_diff` = `beta[hi]` - `beta[lo]`)
    }
  }
    
  return(diff_draws)
}

get_finger <- function(resp, exp_type = NA) {
  if (exp_type == 'norep') {
    finger <- ifelse(resp %in% c('w', 'o'), 'index',
          ifelse(resp %in% c('q', 'p'), 'middle', NA)) 
  } else {
    finger <- ifelse(resp %in% c('left', 'down'), 'index',
          ifelse(resp %in% c('up', 'right'), 'middle', NA))
  }
  return(finger)
}

get_hand <- function(resp, exp_type = NA) {
  if (exp_type == 'norep') {
    hand <- ifelse(resp %in% c('o', 'p'), 'right',
          ifelse(resp %in% c('q', 'w'), 'left', NA))
  } else {
    hand <- ifelse(resp %in% c('left', 'right'), 'right',
        ifelse(resp %in% c('up', 'down'), 'left', NA))
  }
  return(hand)
}

get_lambda <- function(draws) {
  # mixing parameter for rpm_mix
  print("deprecated")
  lambda_draws <- draws %>% 
    subset_draws(variable = 'delta_lambda*', regex = TRUE) %>%
    mutate(lambda = inv_logit(`delta_lambda[1]`))
  return(lambda_draws)
}

get_mu <- function(draws, R_mu = 2, covar = TRUE, transform = inv_logit) {
  raw_draws <- draws 

  if (!covar) {
    if (R_mu == 1) {
      mu_draws <- raw_draws %>% mutate_variables(mu = transform(mu_loc))
    } else {
      mu_draws <- raw_draws %>%
        mutate_variables(
          `mu[h]` = transform(`mu_loc[1]`),
          `mu[g]` = transform(`mu_loc[2]`))
    }
  } else {
    if (R_mu == 1) {
      mu_draws <- raw_draws %>%
        mutate_variables(
          `mu[lo]` = transform(mu_loc + delta_mu_loc * -.5),
          `mu[hi]` = transform(mu_loc + delta_mu_loc * .5))
    } else {
      mu_draws <- raw_draws %>% 
        mutate_variables(
          `mu[h,lo]` = transform(`mu_loc[1]` + `delta_mu_loc[1]` * -.5),
          `mu[h,hi]` = transform(`mu_loc[1]` + `delta_mu_loc[1]` * .5),
          `mu[g,lo]` = transform(`mu_loc[2]` + `delta_mu_loc[2]` * -.5),
          `mu[g,hi]` = transform(`mu_loc[2]` + `delta_mu_loc[2]` * .5))
    }
  }
  return(mu_draws)
}

get_pos_str <- function(pos) {
  pos_str <- 
    ifelse(pos == "[\"0\",\"0.25\"]", 'up', 
      ifelse(pos == "[\"0\",\"-0.25\"]", 'down',
        ifelse(pos == "[\"-0.25\",\"0\"]", 'left',
          ifelse(pos == "[\"0.25\",\"0\"]", 'right', NA))))
  return(pos_str)
}

get_prep1 <- function(draws, transform = inv_logit, time = seq(0, 1, .001)) {

  mus <- get_mu(draws, R_mu = 1, covar = FALSE, transform = transform)
  sigmas <- get_sigma(draws, R_sigma = 1, covar = FALSE, transform = transform)

  prep <- mus %>% 
    inner_join(sigmas) %>% 
    group_by(.draw) %>%
    summarise(dens = dnorm(time, mu, sigma)) %>%
    mutate(t = time) %>%
    group_by(t) %>%
    median_qi()

  return(prep)
}
get_prep2 <- function(draws, transform = inv_logit, time = seq(0, 1, .001)) {
  mus <- get_mu(draws, R_mu = 2, covar = FALSE, transform = transform) %>% spread_draws(mu[r])
  sigmas <- get_sigma(draws, R_sigma = 2, covar = FALSE, transform = transform) %>% spread_draws(sigma[r])

  prep <- mus %>% 
    inner_join(sigmas) %>% 
    group_by(.draw, r) %>%
    summarise(dens = dnorm(time, mu, sigma)) %>%
    mutate(t = time) %>%
    group_by(t, r) %>%
    median_qi()

  return(prep)
}

get_resp_type <- function(resp, goal, exp_type) {
  if (resp == goal) {
    resp_type = 'correct'
  } else {
    if (get_hand(resp, exp_type) == get_hand(goal, exp_type)) {
      resp_type = 'finger_error'
    } else {
      if (get_finger(resp, exp_type) == get_finger(goal, exp_type)) {
        resp_type = 'hand_error'
      } else {
        resp_type = 'full_error'
      }
    }
  }
  return(resp_type)
}

get_sigma <- function(draws, R_sigma = 2, covar = TRUE, transform = inv_logit) {
  # variability in preparation demands
  raw_draws <- draws #%>% subset_draws(variable = 'delta_sigma*', regex = TRUE)

  if (!covar) {
    if (R_sigma == 1) {
      sigma_draws <- raw_draws %>% mutate_variables(sigma = transform(sigma_loc))
    } else {
      sigma_draws <- raw_draws %>%
        mutate_variables(
          `sigma[h]` = transform(`sigma_loc[1]`),
          `sigma[g]` = transform(`sigma_loc[2]`))
    }
  } else {
    if (R_sigma == 1) {
      sigma_draws <- raw_draws %>%
        mutate_variables(
          `sigma[lo]` = transform(sigma_loc + delta_sigma_loc * -.5),
          `sigma[hi]` = transform(sigma_loc + delta_sigma_loc * .5))
    } else {
      sigma_draws <- raw_draws %>% 
        mutate_variables(
          `sigma[h,lo]` = transform(`sigma_loc[1]` + `delta_sigma_loc[1]` * -.5),
          `sigma[h,hi]` = transform(`sigma_loc[1]` + `delta_sigma_loc[1]` * .5),
          `sigma[g,lo]` = transform(`sigma_loc[2]` + `delta_sigma_loc[2]` * -.5),
          `sigma[g,hi]` = transform(`sigma_loc[2]` + `delta_sigma_loc[2]` * .5))
    }
  }
  return(sigma_draws)
}

inv_logit <- function(x) {
  y <- exp(x) / (1 + exp(x))
  return(y)
}
inv_logit1.5 <- function(x) {
  y <- exp(x) / (1 + exp(x))
  return(1.5 * y)
}
inv_logit2 <- function(x) {
  y <- exp(x) / (1 + exp(x))
  return(2 * y)
}

plot_prep1 <- function(draws, mod_name, transform = inv_logit, time = seq(0, 1, .001)) {

  mus <- get_mu(draws, R_mu = 1, covar = FALSE, transform = transform)
  sigmas <- get_sigma(draws, R_sigma = 1, covar = FALSE, transform = transform)

  prep <- mus %>% 
    inner_join(sigmas) %>% 
    group_by(.draw) %>%
    summarise(dens = dnorm(time, mu, sigma)) %>%
    mutate(t = time) %>%
    group_by(t) %>%
    median_qi()

  p <- prep %>% 
    ggplot(aes(x = t*1000, y = dens, ymax = .upper, ymin = .lower)) + 
    geom_ribbon(colour = NA, alpha = .1) +
    geom_line() +
    my_theme +  
    theme(axis.line.y = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())

  path <- sprintf('fig/%s_prep_dens1.jpg', mod_name)
  ggsave(plot = p, filename = path, units = 'in', height = 1.5, width = 3, dpi = 1000)

  return(p)
}

plot_prep1_cov <- function(draws, mod_name, transform = inv_logit, time = seq(0, 1, .001)) {

  mus <- get_mu(draws, R_mu = 1, covar = TRUE, transform = transform) %>% spread_draws(mu[v])
  sigmas <- get_sigma(draws, R_sigma = 1, covar = TRUE, transform = transform) %>% spread_draws(sigma[v])

  prep <- mus %>% 
    inner_join(sigmas) %>% 
    group_by(.draw, v) %>%
    summarise(dens = dnorm(time, mu, sigma)) %>%
    mutate(t = time) %>%
    group_by(t, v) %>%
    median_qi()

  p <- prep %>% 
    ggplot(aes(x = t*1000, y = dens, ymax = .upper, ymin = .lower, fill = v, colour = v, group = v)) + 
    geom_line() +
    geom_ribbon(colour = NA, alpha = .1) +
    scale_fill_manual('cov', values = WONG[c(7, 4)]) +
    scale_colour_manual('cov', values = WONG[c(7, 4)]) +
    my_theme +  
    theme(axis.line.y = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())

  return(p)
}

plot_prep2 <- function(draws, mod_name, transform = inv_logit, time = seq(0, 1, .001)) {
  mus <- draws %>% spread_draws(mu[r])
  sigmas <- draws %>% spread_draws(sigma[r])

  prep <- mus %>% 
    inner_join(sigmas) %>% 
    group_by(.draw, r) %>%
    summarise(dens = dnorm(time, mu, sigma)) %>%
    mutate(t = time, r = as.character(r)) %>%
    group_by(t, r) %>%
    median_qi()

  p <- prep %>%
    ggplot(aes(x = t*1000, y = dens, ymax = .upper, ymin = .lower, fill = r, colour = r, group = r)) + 
    geom_line() +
    geom_ribbon(colour = NA, alpha = .1) +
    scale_fill_manual('Response', values = WONG[c(4, 2)]) +
    scale_colour_manual('Response', values = WONG[c(4, 2)]) +
    my_theme +  
    theme(axis.line.y = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())

  return(p)
}

plot_prep2_cov <- function(draws, mod_name, transform = inv_logit, time = seq(0, 1, .001)) {
  mus <- get_mu(draws, R_mu = 2, covar = TRUE, transform = transform) %>% spread_draws(mu[r, v])
  sigmas <- get_sigma(draws, R_sigma = 2, covar = TRUE, transform = transform) %>% spread_draws(sigma[r, v])

  prep <- mus %>% 
    inner_join(sigmas) %>% 
    group_by(.draw, r, v) %>%
    summarise(dens = dnorm(time, mu, sigma)) %>%
    mutate(t = time) %>%
    group_by(t, r, v) %>%
    median_qi()

  p <- prep %>% mutate(v = factor(v, c('lo', 'hi'))) %>%
    ggplot(aes(x = t*1000, y = dens, ymax = .upper, ymin = .lower, alpha = v, fill = r, colour = r, group = interaction(r, v))) + 
    geom_line() +
    geom_ribbon(colour = NA, alpha = .1) +
    scale_fill_manual('Response', values = WONG[c(4, 2)]) +
    scale_colour_manual('Response', values = WONG[c(4, 2)]) +
    scale_alpha_manual('cov', values = c(.5, 1)) +
    my_theme +  
    theme(axis.line.y = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())

  return(p)
}



se <- function(x, na.rm = TRUE) {
  sd(x, na.rm = TRUE) / sqrt(length(x))
}

smooth_obs <- function(dat, lb = 0, ub = 1000, width = 100, covar_names = c()) {
  # given dataset and model for a single experiment
  # run sliding window analysis on accuracy over time
  # return smoothed time-series of accuracy
  obs <- dat
  tmp <- obs

  if (length(covar_names) > 0) {
    obs$covar1 <- obs[[covar_names[1]]]
    tmp <- obs %>% group_by(covar1)
  }
  if (length(covar_names) > 1) {
    obs$covar2 <- obs[[covar_names[2]]]
    tmp <- obs %>% group_by(covar1, covar2)
  }
  if (length(covar_names) > 2) {
    obs$covar3 <- obs[[covar_names[3]]]
    tmp <- obs %>% group_by(covar1, covar2, covar3)
  }

  sw_dat <- tmp %>% summarise(
    PT = lb:ub, 
    obs_mu = sw_smooth(x = X, y = Y, width = width, lb = lb, ub = ub, FUN = mean),
    obs_se = sw_smooth(x = X, y = Y, width = width, lb = lb, ub = ub, FUN = se))

  return(sw_dat)
}

smooth_obs_pred <- function(dat, fit, lb = 0, ub = 1000, width = 100, covar_names = c()) {
  # given dataset and model for a single experiment
  # run sliding window analysis on accuracy over time
  # return smoothed time-series of accuracy
  
  pred <- t(rstan::extract(fit, pars = 'theta')$theta)

  obs_pred <- dat %>% 
    mutate(
      theta = apply(X = pred, FUN = mean, MARGIN = 1),
      high = apply(X = pred, FUN = ci95, MARGIN = 1),
      low = apply(X = pred, FUN = ci05, MARGIN = 1))

  if (length(covar_names) > 0) {
    obs_pred$covar1 <- obs_pred[[covar_names[1]]]
    tmp <- obs_pred %>% group_by(covar1)
  }
  if (length(covar_names) > 1) {
    obs_pred$covar2 <- obs_pred[[covar_names[2]]]
    tmp <- obs_pred %>% group_by(covar1, covar2)
  }
  if (length(covar_names) > 2) {
    obs_pred$covar3 <- obs_pred[[covar_names[3]]]
    tmp <- obs_pred %>% group_by(covar1, covar2, covar3)
  }

  sw_dat <- tmp %>% summarise(
    PT = lb:ub, 
    obs_mu = sw_smooth(x = X, y = Y, lb = lb, ub = ub),
    obs_se = sw_smooth(x = X, y = Y, lb = lb, ub = ub, FUN = se),
    pred  = sw_smooth(x = X, y = theta, lb = lb, ub = ub),
    upper = sw_smooth(x = X, y = high, lb = lb, ub = ub),
    lower = sw_smooth(x = X, y = low, lb = lb, ub = ub))

  return(sw_dat)
}

smooth_pred <- function(dat, lb = 0, ub = 1000, width = 100) {
  # given dataset and model for a single experiment
  # run sliding window analysis on accuracy over time
  # return smoothed time-series of accuracy

  sw_dat <- dat %>% summarise(
    PT = lb:ub, 
    pred  = sw_smooth(x = X, y = theta, lb = lb, ub = ub),
    upper = sw_smooth(x = X, y = high, lb = lb, ub = ub),
    lower = sw_smooth(x = X, y = low, lb = lb, ub = ub))

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

    if (!is_empty(sw)) {
      sw <- append(sw, FUN(win, na.rm = TRUE))
    } else {
      sw <- append(sw, NA) 
    }

  }
  return(sw)
}

sw_ttest <- function(dat, lb = 0, ub = 1000, width = 100, covar_name) {
  # ** WARNING: VERY LONG RUNTIME! **

  # perform multilevel bernoulli regression with binary predictor over sliding window

  my_priors <- c(set_prior("normal(0, 1)", class = "b"))

  sw <- c()

  for (i in lb:ub) {

    lower <- i - (width / 2) 
    upper <- i + (width / 2)

    win <- dat %>% filter(PT <= upper & PT >= lower)

    if (!is_empty(sw)) {

      if (i == lb) {

        if (covar_name == 'rew') {
          m_formula <- bf(y ~ rew + (rew | id))
        } else {
          if (covar_name == 'congr') {
            m_formula <- bf(y ~ congr + (congr | id))
          } else {
            if (covar_name == 'congr0') {
              m_formula <- bf(y ~ congr0 + (congr0 | id))
            } else {
              if (covar_name == 'err0') {
                m_formula <- bf(y ~ err0 + (err0 | id))
              }
            }
          }
        }
        
        m <- brm(formula = m_formula, data = win, prior = my_priors, family = bernoulli, cores = 1)

      } else {

        m <- update(m, newdata = win,  cores = 1)

      }

      p_dir <- m %>% p_direction() %>% tibble() %>% filter(Parameter == "b_rew") %>% select(pd)
      eff <- m %>% spread_draws(b_rew) %>% median_qi() %>% select(b_rew)
      dp_dir <- p_dir * sign(eff)
      sw <- append(sw, dp_dir$pd)

    } else {

      sw <- append(sw, NA) 

    }

  }

  return(sw)

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
    axis.line = element_line(size = .5)))


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
      values = c(WONG[2], WONG[3], WONG[7], WONG[6])),
    scale_fill_manual(
      values = c(WONG[2], WONG[3], WONG[7], WONG[6])),
    scale_alpha_manual(
      values = c(.5, 1)),
    coord_cartesian(
      xlim = c(0, 1000), 
      ylim = c(0, 1)))

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
