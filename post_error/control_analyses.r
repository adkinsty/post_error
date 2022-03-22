library(tidyverse)
library(tidybayes)


covar_name <- "err0"
mod_type <- "rpm1_hier"

general_path <- "/Users/adkinsty/Dropbox (University of Michigan)/"
specific_path <- "LeeLab/Experiments/Exp_files/forced_response/post_error/"
setwd(sprintf("%s%s", general_path, specific_path))
source("../../forced_response/func.R")

err0_model <- "rpm1_hier_mu_err0_sigma_err0_beta_err0"
err00_model <- "rpm1_hier_mu_err00_sigma_err00_beta_err00"

exp1 <- "memory"
exp2 <- "memory2"
exp3 <- "memory_color_short"
exp4 <- "memory_color_long"

dat1 <- read_csv(sprintf("%s/dat.csv", exp1)) %>%
    select(PT, err0, err00, accuracy, key_rep, key_stim_rep, participant) %>%
    mutate(exp = exp1)
dat2 <- read_csv(sprintf("%s/dat.csv", exp2)) %>%
    select(PT, err0, err00, accuracy, key_rep, key_stim_rep, participant) %>%
    mutate(exp = exp2)
dat3 <- read_csv(sprintf("%s/dat.csv", exp3)) %>%
    select(PT, err0, err00, accuracy, key_rep, key_stim_rep, participant) %>%
    mutate(exp = exp3)
dat4 <- read_csv(sprintf("%s/dat.csv", exp4)) %>%
    select(PT, err0, err00, accuracy, key_rep, key_stim_rep, participant) %>%
    mutate(exp = exp4)


 ## Stay/shift after error
 rep_summary <- dat1  %>%
    rbind(dat2) %>%
    rbind(dat3) %>%
    rbind(dat4) %>%
    filter(PT > 500) %>%
    filter(accuracy == 0) %>%
    group_by(err0, exp) %>%
    summarise(
        rep_mean = mean(key_rep),
        rep_upper = rep_mean + se(key_rep) * 1.96,
        rep_lower = rep_mean - se(key_rep) * 1.96,
        n = n()
)

dodge_pos <- position_dodge(.25)
rep_summary %>%
    mutate(
        exp = factor(exp, levels = c(
            "memory", "memory2", "memory_color_short", "memory_color_long"
))) %>%
    ggplot(aes(
        x = as.character(err0),
        y = rep_mean,
        ymax = rep_upper,
        ymin = rep_lower,
        colour = exp, group = exp)) +
    geom_point(position = dodge_pos) +
    geom_line(position = dodge_pos) +
    geom_errorbar(width = 0, position = dodge_pos) +
    scale_y_continuous("P(repeat)") +
    scale_x_discrete("N-1 Accuracy", labels = c("correct", "error")) +
    scale_colour_manual("Experiment",
        labels = 1:4,
        values = WONG[c(5, 8, 3, 2)]) +
    theme_classic(base_size = 10) +
    coord_cartesian(ylim = c(0, 1))
ggsave("fig/key_rep.jpg", dpi = 1000, units = "in", height = 2, width = 3)


priors <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior("normal(0, 1)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "sd"))

fit_rep1 <- brm(
  data = dat1 %>% filter(PT > 500) %>% filter(accuracy == 0),
  formula = key_rep ~ err0 + (err0 | participant),
  family = "bernoulli",
  prior = priors,
  chains = 4, cores = 1, iter = 4000, warmup = 2000)

fit_rep2 <- update(fit_rep1,
    newdata = dat2 %>% filter(PT > 500) %>% filter(accuracy == 0))
fit_rep3 <- update(fit_rep1,
    newdata = dat3 %>% filter(PT > 500) %>% filter(accuracy == 0))
fit_rep4 <- update(fit_rep1,
    newdata = dat4 %>% filter(PT > 500) %>% filter(accuracy == 0))





## ITI effect
fit3 <- read_rds(sprintf("%s/fit/%s_%s.rds", exp3, exp3, err0_model))
draws3 <- as_draws_df(fit3)
betas3 <- get_beta(draws3, R_beta = 1, covar = TRUE)
beta_diff3 <- get_diff(betas3, par = "beta", R = 1) %>%
    subset_draws(variable = "beta_diff")
beta_diff3 %>% median_qi()

fit4 <- read_rds(sprintf("%s/fit/%s_%s.rds", exp4, exp4, err0_model))
draws4 <- as_draws_df(fit4)
betas4 <- get_beta(draws4, R_beta = 1, covar = TRUE)
beta_diff4 <- get_diff(betas4, par = "beta", R = 1) %>%
    subset_draws(variable = "beta_diff")
beta_diff4 %>% median_qi()

iti_dat <- tibble(
    diff = c(beta_diff3$beta_diff, beta_diff4$beta_diff),
    exp = rep(c("short", "long"), each = nrow(beta_diff3)))

iti_dat %>%
    ggplot(aes(x = diff, group = exp, fill = exp, colour = exp)) +
    geom_histogram(
        binwidth = .005,
        position = "identity",
        alpha = .5) +
    scale_x_continuous(expression(Delta[beta])) +
    scale_y_continuous("MCMC Draws") +
    scale_fill_manual("ITI", values = WONG[c(2, 3)]) +
    scale_colour_manual("ITI", values = WONG[c(2, 3)]) +
    theme_classic(base_size = 10)
ggsave("fig/iti.jpg", dpi = 1000, units = "in", height = 2, width = 3)




## N-2 err effect after N-1 correct
dat1 %>%
    filter(!is.na(err0) & !is.na(err00)) %>%
    mutate(X = PT, Y = accuracy) %>%
    smooth_obs(covar_names = c("err00", "err0"), width = 100) %>%
    rename(err00 = covar1, err0 = covar2) %>%
    filter(err0 < 0) %>%
    ggplot(aes(
        x = PT,
        y = obs_mu,
        ymax = obs_mu + obs_se * 1.96,
        ymin = obs_mu - obs_se * 1.96,
        colour = as.character(err00),
        fill = as.character(err00))) +
    sw_geom +
    scale_x_continuous("PT (ms)",
        breaks = seq(0, 1000, 100), labels = as.character(seq(0, 2000, 200))) +
    scale_colour_manual(values = WONG[c(4, 7)]) +
    scale_fill_manual(values = WONG[c(4, 7)]) +
    my_theme
ggsave("fig/err00.jpg", dpi = 1000, units = "in", height = 2, width = 3.5)
 

stan_fit <- read_rds(sprintf("%s/fit/%s_%s.rds", exp1, exp1, err00_model))
stan_fit <- read_rds(sprintf("%s/fit/%s_%s.rds", exp2, exp2, err00_model))
stan_fit <- read_rds(sprintf("%s/fit/%s_%s.rds", exp3, exp3, err00_model))
stan_fit <- read_rds(sprintf("%s/fit/%s_%s.rds", exp4, exp4, err00_model))

stan_fit <- read_rds(sprintf("%s/fit/%s_err0_%s.rds", exp1, exp1, err00_model))
stan_fit <- read_rds(sprintf("%s/fit/%s_err0_%s.rds", exp2, exp2, err00_model))
stan_fit <- read_rds(sprintf("%s/fit/%s_err0_%s.rds", exp3, exp3, err00_model))
stan_fit <- read_rds(sprintf("%s/fit/%s_err0_%s.rds", exp4, exp4, err00_model))


draws <- as_draws_df(stan_fit)
betas <- get_beta(draws, R_beta = 1, covar = TRUE)
beta_diff <- get_diff(betas, par = "beta", R = 1) %>%
    subset_draws(variable = "beta_diff") %>%
    median_qi()
