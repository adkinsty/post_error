// response preparation model with flexible inputs
// inputs determine number of processes and covariates for each param
functions {
  vector get_theta(int N, vector time, int R_mu, int R_sigma, matrix x_mu, matrix x_sigma, matrix x_beta, matrix x_alpha, matrix x_lambda, matrix delta_mu,  matrix delta_sigma, vector delta_beta, vector delta_alpha, vector delta_lambda) {
    
    matrix[N, R_mu] mu;
    matrix[N, R_sigma] sigma;
    vector[N] beta;
    vector[N] alpha;
    vector[N] lambda;

    vector[N] theta; 
    vector[1] phi; 
    vector[1] neg_phi;
    vector[2] psi; 

    mu = inv_logit(x_mu * delta_mu);
    sigma = inv_logit(x_sigma * delta_sigma);
    beta = inv_logit(x_beta * delta_beta);
    alpha = inv_logit(x_alpha * delta_alpha);
    lambda = inv_logit(x_lambda * delta_lambda);

    for (n in 1:N) {

      phi[1] = lambda[n] * normal_cdf(time[n], mu[n, 1], sigma[n, 1])  +  (1 - lambda[n]) * normal_cdf(time[n], mu[n, 2], sigma[n, 2]);
      neg_phi[1] = 1 - phi[1];

      psi[1] = neg_phi[1] * alpha[n]; 
      psi[2] = phi[1]  * beta[n];

      theta[n] = sum(psi);  // response prob (@t)
    }
    return(theta);
  }
  
}
data {
  
  int<lower=1> N; // n observations

  int<lower=0,upper=1> y[N]; // observations
  vector<lower=0,upper=1>[N] time; // time in seconds

  int<lower=1> R_mu; // n process for mu
  int<lower=1> R_sigma; // n process for sigma

  int<lower=1> K_mu; // n covariates for mu
  int<lower=1> K_sigma; // n covariates for sigma
  int<lower=1> K_beta; // n covariates for beta
  int<lower=1> K_alpha; // n covariates for alpha
  int<lower=1> K_lambda; // n covariates for lambda

  matrix[N, K_mu] x_mu; // design matrix for mu
  matrix[N, K_sigma] x_sigma; // design matrix for sigma
  matrix[N, K_beta] x_beta; // design matrix for beta
  matrix[N, K_alpha] x_alpha; // design matrix for alpha
  matrix[N, K_alpha] x_lambda; // design matrix for lambda

  int<lower=0,upper=1> gq; 
  int<lower=0,upper=1> prior_only;

}
parameters {

  matrix[K_mu, R_mu] delta_mu;
  matrix[K_sigma, R_sigma] delta_sigma;
  vector[K_beta] delta_beta;
  vector[K_alpha] delta_alpha;
  vector[K_lambda] delta_lambda;

}
transformed parameters {
  vector[N] theta;
  theta = get_theta(N, time, R_mu, R_sigma, x_mu, x_sigma, x_beta, x_alpha, x_lambda, delta_mu, delta_sigma, delta_beta, delta_alpha, delta_lambda);
}
model {

  // priors for intercepts
  delta_mu[1, 1] ~ normal(-1, .5);
  delta_sigma[1, 1] ~ normal(-2, .5);
  delta_beta[1] ~ normal(1, .5);
  delta_alpha[1] ~ normal(-1, .5);

  // priors for additional process intercepts, if applicable
  if (R_mu > 1) {
    delta_mu[1, R_mu] ~ normal(0, .5);
  }
  if (R_sigma > 1) {
    delta_sigma[1, R_sigma] ~ normal(-2, .5);
  }

  // prior for mixing parameter
  delta_lambda[1] ~ normal(-2, .5);

  // priors for slopes, if applicable
  if (K_mu > 1) {
    to_vector(delta_mu[2:K_mu, ]) ~ normal(0, .5);
  }
  if (K_sigma > 1) {
    to_vector(delta_sigma[2:K_sigma, ]) ~ normal(0, .5);
  }
  if (K_beta > 1) {
    delta_beta[2:K_beta] ~ normal(0, .5);
  }
  if (K_alpha > 1) {
    delta_alpha[2:K_alpha] ~ normal(0, .5);
  }
  if (K_lambda > 1) {
    delta_lambda[2:K_lambda] ~ normal(0, .5);
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
