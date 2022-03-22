// response preparation model with flexible inputs
// inputs determine number of processes and covariates for each param
functions {
  vector get_theta(int h, int g, int N, vector time, vector congr, int R_mu, int R_sigma, int R_beta, matrix x_mu, matrix x_sigma, matrix x_beta, matrix x_alpha, matrix delta_mu,  matrix delta_sigma, matrix delta_beta, vector delta_alpha) {
    
    matrix[N, R_mu] mu;  // h prep mean, g prep mean
    matrix[N, R_sigma] sigma; // h prep sd, g prep sd
    matrix[N, R_beta] beta; // p(corr | h prep), p(corr | g prep)
    vector[N] alpha; // p(corr | not h prep & not g prep)

    vector[2] phi; // p(h prep | t),  p(g prep | t)
    vector[2] not_phi; // // p(not h prep),  p(not g prep)
    vector[4] psi; // p(corr | not h prep & not g prep),  p(corr | h prep & not g prep), p(corr | not h prep & g prep), p(corr | h prep & g prep)
  
    vector[N] theta; // p(corr | t)

    mu = inv_logit(x_mu * delta_mu);
    sigma = inv_logit(x_sigma * delta_sigma);
    beta = inv_logit(x_beta * delta_beta);
    alpha = inv_logit(x_alpha * delta_alpha);

    for (n in 1:N) {

      phi[h] = normal_cdf(time[n], mu[n, h], sigma[n, h]); // prob. that habit is prepared
      phi[g] = normal_cdf(time[n], mu[n, R_mu], sigma[n, R_mu]); // prob. that goal is prepared
      not_phi[h] = 1 - phi[h];
      not_phi[g] = 1 - phi[g];

      psi[1] = not_phi[h] * not_phi[g] * alpha[n]; 
      
      if (congr[n] > 0) {
        psi[2] = phi[h]  * not_phi[g] * beta[n, h]; 
      } else {
        psi[2] = phi[h]  * not_phi[g] * (1 - beta[n, h]); 
      }
      
      psi[3] = not_phi[h] * phi[g] * beta[n, R_beta];
      psi[4] = phi[h]  * phi[g] * beta[n, R_beta];

      theta[n] = sum(psi);  // response prob (@t)
    }
    return(theta);
  }
  
}
data {
  
  int<lower=1> N; // n observations

  int<lower=0, upper=1> y[N]; // observations
  vector<lower=0, upper=1>[N] time; // time in seconds
  vector[N] congr; // S-R congruency switch

  int<lower=1, upper=2> R_mu; // n process for mu
  int<lower=1, upper=2> R_sigma; // n process for sigma
  int<lower=1, upper=2> R_beta; // n process for beta

  int<lower=1> K_mu; // n covariates for mu
  int<lower=1> K_sigma; // n covariates for sigma
  int<lower=1> K_beta; // n covariates for beta
  int<lower=1> K_alpha; // n covariates for alpha

  matrix[N, K_mu] x_mu; // design matrix for mu
  matrix[N, K_sigma] x_sigma; // design matrix for sigma
  matrix[N, K_beta] x_beta; // design matrix for beta
  matrix[N, K_alpha] x_alpha; // design matrix for alpha

  int<lower=0,upper=1> gq; 
  int<lower=0,upper=1> prior_only;

}
transformed data {
  int h;
  int g;

  h = 1;
  g = 2;
}
parameters {

  matrix[K_mu, R_mu] delta_mu;
  matrix[K_sigma, R_sigma] delta_sigma;
  matrix[K_beta, R_beta] delta_beta;
  vector[K_alpha] delta_alpha;

}
transformed parameters {
  vector[N] theta;
  theta = get_theta(h, g, N, time, congr, R_mu, R_sigma, R_beta, x_mu, x_sigma, x_beta, x_alpha, delta_mu, delta_sigma, delta_beta, delta_alpha);
}
model {

  // priors for intercepts
  delta_mu[1, h] ~ normal(-1, .5);
  delta_sigma[1, h] ~ normal(-1, .5);
  delta_beta[1, h] ~ normal(1, .5);
  delta_alpha[1] ~ normal(0, .5);
  
  // priors for additional process intercepts, if applicable
  if (R_mu > 1) {
    delta_mu[1, g] ~ normal(0, .5); // bias mu_g to be greater than mu_h
  }
  if (R_sigma > 1) {
    delta_sigma[1, g] ~ normal(-1, .5);
  }
  if (R_beta > 1) {
    delta_beta[1, g] ~ normal(1, .5);
  }

  // priors for slopes, if applicable
  if (K_mu > 1) {
    to_vector(delta_mu[2:K_mu, ]) ~ normal(0, .5);
  }
  if (K_sigma > 1) {
    to_vector(delta_sigma[2:K_sigma, ]) ~ normal(0, .5);
  }
  if (K_beta > 1) {
    to_vector(delta_beta[2:K_beta, ]) ~ normal(0, .5);
  }
  if (K_alpha > 1) {
    delta_alpha[2:K_alpha] ~ normal(0, .5);
  }

  if (prior_only == 0) {
    y ~ bernoulli(theta);
  }

}
generated quantities {

  vector[N * gq] log_lik;

  if (gq == 1) {

    for (n in 1:N) {
      
      log_lik[n] = bernoulli_lpmf(y[n] | theta[n]);

    }
  }
}
