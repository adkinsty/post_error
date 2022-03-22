// response preparation model with flexible inputs
// inputs determine number of processes and covariates for each param
functions {
  vector get_theta(int N, vector time, vector congr, int R_mu, int R_sigma, int R_beta, matrix x_mu, matrix x_sigma, matrix x_beta, matrix x_alpha, matrix delta_mu,  matrix delta_sigma, matrix delta_beta, vector delta_alpha) {
    
    matrix[N, R_mu] mu;
    matrix[N, R_sigma] sigma;
    matrix[N, R_beta] beta;
    vector[N] alpha;
    vector[N] theta; 
    vector[2] phi; 
    vector[2] neg_phi;
    vector[4] psi; 

    mu = inv_logit(x_mu * delta_mu);
    sigma = inv_logit(x_sigma * delta_sigma);
    beta = inv_logit(x_beta * delta_beta);
    alpha = inv_logit(x_alpha * delta_alpha);

    for (n in 1:N) {

      phi[1] = normal_cdf(time[n], mu[n, 1], sigma[n, 1]); // prob. that habit is prepared
      phi[2] = normal_cdf(time[n], mu[n, R_mu], sigma[n, R_sigma]); // prob. that goal is prepared
      neg_phi[1] = 1 - phi[1];
      neg_phi[2] = 1 - phi[2];

      psi[1] = neg_phi[1] * neg_phi[2] * alpha[n]; 
      
      if (congr[n] > 0) {
        psi[2] = phi[1]  * neg_phi[2] * beta[n, 1]; 
      } else {
        psi[2] = phi[1]  * neg_phi[2] * (1 - beta[n, 1]); 
      }
      
      psi[3] = neg_phi[1] * phi[2] * beta[n, R_beta];
      psi[4] = phi[1]  * phi[2] * beta[n, R_beta];

      theta[n] = sum(psi);  // response prob (@t)
    }
    return(theta);
  }
  
}
data {
  
  int<lower=1> N; // n observations

  int<lower=0,upper=1> y[N]; // observations
  vector<lower=0,upper=1>[N] time; // time in seconds
  vector[N] congr; // S-R congruency switch

  int<lower=1> R_mu; // n process for mu
  int<lower=1> R_sigma; // n process for sigma
  int<lower=1> R_beta; // n process for beta

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

  matrix[K_mu, R_mu] delta_mu;
  matrix[K_sigma, R_sigma] delta_sigma;
  matrix[K_beta, R_beta] delta_beta;
  vector[K_alpha] delta_alpha;

}
transformed parameters {
  vector[N] theta;
  theta = get_theta(N, time, congr, R_mu, R_sigma, R_beta, x_mu, x_sigma, x_beta, x_alpha, delta_mu, delta_sigma, delta_beta, delta_alpha);
}
model {

  // informative priors for intercepts, to improve identifiability
  delta_mu[1, 1] ~ normal(-2, .25);
  delta_sigma[1, 1] ~ normal(-2, .25);
  delta_beta[1, 1] ~ normal(1, .25);
  delta_alpha[1] ~ normal(0, .25);
  
  // priors for additional process's parameters, if applicable
  if (R_mu > 1) {
    delta_mu[1, R_mu] ~ normal(-1, .25);
  }
  if (R_sigma > 1) {
    delta_sigma[1, R_sigma] ~ normal(-1, .25);
  }
  if (R_beta > 1) {
    delta_beta[1, R_beta] ~ normal(2, .25);
  }

  // weakly informative priors for covariate effects, if applicable
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
