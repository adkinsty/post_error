// response preparation model with flexible inputs
// inputs determine number of processes and covariates for each param
functions {
  vector get_theta(int N, vector time, matrix x_mu, matrix x_sigma, matrix x_beta, matrix x_alpha, vector delta_mu,  vector delta_sigma, vector delta_beta, vector delta_alpha) {
    
    vector[N] mu; // prep mean
    vector[N] sigma; // prep sd
    vector[N] beta; // p(corr | prep)
    vector[N] alpha; // p(corr | not prep)

    real phi; // p(prep | t)
    real not_phi; // p(not prep | t)
    vector[2] psi;  // p(corr | prep) & p(corr | not prep) 

    vector[N] theta; // p(corr | t)

    mu = inv_logit(x_mu * delta_mu);
    sigma = inv_logit(x_sigma * delta_sigma);
    beta = inv_logit(x_beta * delta_beta);
    alpha = inv_logit(x_alpha * delta_alpha);

    for (n in 1:N) {

      phi = normal_cdf(time[n], mu[n], sigma[n]);
      not_phi = 1 - phi;

      psi[1] = not_phi * alpha[n]; 
      psi[2] = phi * beta[n];

      theta[n] = sum(psi);  // response prob (@t)
    }
    return(theta);
  }
  
}
data {
  
  int<lower=1> N; // n observations

  int<lower=0,upper=1> y[N]; // observations
  vector<lower=0,upper=1>[N] time; // time in seconds

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
parameters {

  vector[K_mu] delta_mu;
  vector[K_sigma] delta_sigma;
  vector[K_beta] delta_beta;
  vector[K_alpha] delta_alpha;

}
transformed parameters {
  vector[N] theta;
  theta = get_theta(N, time, x_mu, x_sigma, x_beta, x_alpha, delta_mu, delta_sigma, delta_beta, delta_alpha);
}
model {

  // priors for intercepts
  delta_mu[1] ~ normal(-1, .5);
  delta_sigma[1] ~ normal(-2, .5);
  delta_beta[1] ~ normal(1, .5);
  delta_alpha[1] ~ normal(-1, .5);

  // priors for slopes, if applicable
  if (K_mu > 1) {
    to_vector(delta_mu[2:K_mu]) ~ normal(0, .5);
  }
  if (K_sigma > 1) {
    to_vector(delta_sigma[2:K_sigma]) ~ normal(0, .5);
  }
  if (K_beta > 1) {
    to_vector(delta_beta[2:K_beta]) ~ normal(0, .5);
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
