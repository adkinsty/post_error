setwd('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/memory/demo')


get_theta <- function(phi, alpha, beta) {
  psi <- c()
  psi[1] = (1 - phi) * alpha
  psi[2] = phi * beta
  return(sum(psi))
}

N <- 1e4
K <- 2

plot_all <- list(geom_line(size = .5), 
  scale_x_continuous("Time", breaks = seq(0, 1, .1)),
  theme_classic(base_size = 10),
  theme(legend.position = "none", axis.text = element_blank()))
plot_dens <- list(plot_all,
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 5)), 
  scale_y_continuous("Density"))
plot_phi <- list(plot_all,
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)), 
  scale_y_continuous(expression(phi), breaks = seq(0, 1, .1)))
plot_theta <- list(plot_all,
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)), 
  scale_y_continuous("Accuracy", breaks = seq(0, 1, .1)))


dat <- tibble(
    alpha = .25, 
    mu0 = .4, 
    sigma0 = .15, 
    beta0 = .9,
    time = runif(N * K), 
    covar = rep(c(-.5, .5), N),
    delta_mu = .05,
    delta_sigma = .05, 
    delta_beta = .05)

dat_mu <- dat %>%
  mutate(mu = mu0 + delta_mu * covar) %>%
  rowwise() %>% 
  mutate(
    dens = dnorm(time, mu, sigma0),
    phi = pnorm(time, mu, sigma0),
    theta = get_theta(phi, alpha, beta0))

dat_mu %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = dens, colour = covar, group = covar)) +
  scale_colour_manual("covar", values = WONG[c(4, 7)], labels = c("A", "B")) +
  plot_dens +
  ggsave("fig/mu_dens.jpg", units = "in", height = 1.5, width = 2, dpi = 1000)

dat_mu %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = phi, colour = covar, group = covar)) +
  scale_colour_manual("covar", values = WONG[c(4, 7)], labels = c("A", "B")) +
  plot_phi +
  ggsave("fig/mu_phi.jpg", units = "in", height = 1.5, width = 2, dpi = 1000)

dat_mu %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = theta, colour = covar, group = covar)) +
  scale_colour_manual("covar", values = WONG[c(4, 7)], labels = c("A", "B")) +
  plot_theta +
  geom_hline(yintercept = .25, linetype = 'dashed', colour = WONG[1], alpha = .5) +
  ggsave("fig/mu_theta.jpg", units = "in", height = 1.5, width = 2, dpi = 1000)




dat_sigma <- dat %>%
  mutate(sigma = sigma0 + delta_sigma * covar) %>%
  rowwise() %>% 
  mutate(
    dens = dnorm(time, mu0, sigma),
    phi = pnorm(time, mu0, sigma),
    theta = get_theta(phi, alpha, beta0))

dat_sigma %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = dens, colour = covar, group = covar)) +
  scale_colour_manual("covar", values = WONG[c(4, 7)], labels = c("A", "B")) +
  plot_dens +
  ggsave("fig/sigma_dens.jpg", units = "in", height = 1.5, width = 2, dpi = 1000)

dat_sigma %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = phi, colour = covar, group = covar)) +
  scale_colour_manual("covar", values = WONG[c(4, 7)], labels = c("A", "B")) +
  plot_phi +
  ggsave("fig/sigma_phi.jpg", units = "in", height = 1.5, width = 2, dpi = 1000)

dat_sigma %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = theta, colour = covar, group = covar)) +
  scale_colour_manual("covar", values = WONG[c(4, 7)], labels = c("A", "B")) +
  plot_theta +
  geom_hline(yintercept = .25, linetype = 'dashed', colour = WONG[1], alpha = .5) +
  ggsave("fig/sigma_theta.jpg", units = "in", height = 1.5, width = 2, dpi = 1000)


dat_beta <- dat %>%
  mutate(beta = beta0 - delta_beta * covar) %>%
  rowwise() %>% 
  mutate(
    dens = dnorm(time, mu0, sigma0),
    phi = pnorm(time, mu0, sigma0),
    theta = get_theta(phi, alpha, beta))

dat_beta %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = dens, colour = covar, group = covar, alpha = .25)) +
  scale_colour_manual("covar", values = WONG[c(4, 7)], labels = c("A", "B")) +
  plot_dens +
  ggsave("fig/beta_dens.jpg", units = "in", height = 1.5, width = 2, dpi = 1000)

dat_beta %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = phi, colour = covar, group = covar, alpha = .25)) +
  scale_colour_manual("covar", values = WONG[c(4, 7)], labels = c("A", "B")) +
  plot_phi +
  ggsave("fig/beta_phi.jpg", units = "in", height = 1.5, width = 2, dpi = 1000)

dat_beta %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = theta, colour = covar, group = covar)) +
  scale_colour_manual("covar", values = WONG[c(4, 7)], labels = c("A", "B")) +
  plot_theta +
  geom_hline(yintercept = .25, linetype = 'dashed', colour = WONG[1], alpha = .5) +
  ggsave("fig/beta_theta.jpg", units = "in", height = 1.5, width = 2, dpi = 1000)



dat_mu_sigma <- dat %>%
  mutate(
    mu = mu0 + delta_mu * covar, 
    sigma = sigma0 + delta_sigma * covar) %>%
  rowwise() %>% 
  mutate(
    dens = dnorm(time, mu, sigma),
    phi = pnorm(time, mu, sigma),
    theta = get_theta(phi, alpha, beta))

dat_mu_sigma %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = dens, colour = covar, group = covar)) +
  scale_colour_manual("covar", values = WONG[c(4, 8)], labels = c("A", "B")) +
  plot_dens +
  ggsave("fig/mu_sigma_dens.jpg", units = "in", height = 3, width = 4, dpi = 1000)

dat_mu_sigma %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = phi, colour = covar, group = covar)) +
  scale_colour_manual("covar", values = WONG[c(4, 8)], labels = c("A", "B")) +
  plot_phi +
  ggsave("fig/mu_sigma_phi.jpg", units = "in", height = 3, width = 4, dpi = 1000)

dat_mu_sigma %>% 
  mutate(covar = as.character(covar)) %>%
  ggplot(aes(x = time, y = theta, colour = covar, group = covar)) +
  scale_colour_manual("covar", values = WONG[c(4, 8)], labels = c("A", "B")) +
  plot_theta +
  ggsave("fig/mu_sigma_theta.jpg", units = "in", height = 3, width = 4, dpi = 1000)

