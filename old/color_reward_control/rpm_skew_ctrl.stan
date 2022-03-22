// response preparation model with flexible inputs
// inputs determine number of processes and covariates for each param
functions {
  vector get_theta(int N, vector time, vector phase, int R_mu, int R_sigma, int R_lambda, int R_beta, matrix x_mu,  matrix x_sigma, matrix x_lambda, matrix x_beta, matrix x_alpha, matrix delta_mu,  matrix delta_sigma, matrix delta_lambda, matrix delta_beta, vector delta_alpha) {
    // transformed parameters defined here to avoid memory overlaod
    matrix[N, R_mu] mu;
    matrix[N, R_sigma] sigma;
    matrix[N, R_lambda] lambda;

    matrix[N, R_mu] xi;
    matrix[N, R_sigma] omega;
    matrix[N, R_lambda] gamma;

    matrix[N, R_beta] beta;
    vector[N] alpha;

    vector[N] theta; 
    real phi; 
    real neg_phi;
    vector[2] psi; 

    mu = inv_logit(x_mu * delta_mu);
    sigma = inv_logit(x_sigma * delta_sigma);
    lambda = x_lambda * delta_lambda;

    // adapted from brms distributions.R
    gamma = lambda ./ sqrt(1 + square(lambda));
    omega = sigma ./ sqrt(1 - to_matrix(rep_array(2 / pi(), N, R_lambda)) .* square(gamma));
    xi = mu - omega .* gamma .* to_matrix(rep_array(sqrt(2 / pi()), N, R_lambda));

    beta = inv_logit(x_beta * delta_beta);
    alpha = inv_logit(x_alpha * delta_alpha);

    for (n in 1:N) {

      if (phase[n] > 0) {
        phi = skew_normal_cdf(time[n], xi[n, 1], omega[n, 1], lambda[n, 1]); // prob. that response to location is prepared
        psi[1] = (1 - phi) * alpha[n]; 
        psi[2] = phi * beta[n, 1]; 
      } else {
        phi = skew_normal_cdf(time[n], xi[n, R_mu], omega[n, R_sigma], lambda[n, R_lambda]); // prob. that response to color is prepared
        psi[1] = (1 - phi) * alpha[n]; 
        psi[2] = phi * beta[n, R_beta]; 
      }
      theta[n] = sum(psi);  // response prob (@t)
    }
    return(theta);
  }
  
}
data {
  
  int<lower=1> N; // n observations

  int<lower=0,upper=1> y[N]; // observations
  vector<lower=0,upper=1>[N] time; // time in seconds 
  vector[N] phase; // location phase (-.5) or color phase (.5)

  int<lower=1> R_mu; // n process for mu
  int<lower=1> R_sigma; // n process for sigma
  int<lower=1> R_lambda; // n process for lambda
  int<lower=1> R_beta; // n process for beta

  int<lower=1> K_mu; // n covariates for mu
  int<lower=1> K_sigma; // n covariates for sigma
  int<lower=1> K_lambda; // n covariates for lambda
  int<lower=1> K_beta; // n covariates for beta
  int<lower=1> K_alpha; // n covariates for alpha

  matrix[N, K_mu] x_mu; // design matrix for mu
  matrix[N, K_sigma] x_sigma; // design matrix for sigma
  matrix[N, K_lambda] x_lambda; // design matrix for lambda
  matrix[N, K_beta] x_beta; // design matrix for beta
  matrix[N, K_alpha] x_alpha; // design matrix for alpha

  int<lower=0,upper=1> gq; 
  int<lower=0,upper=1> prior_only;

}
parameters {

  matrix[K_mu, R_mu] delta_mu;
  matrix[K_sigma, R_sigma] delta_sigma;
  matrix[K_lambda, R_lambda] delta_lambda;
  matrix[K_beta, R_beta] delta_beta;
  vector[K_alpha] delta_alpha;

}
transformed parameters {

  vector[N] theta;

  theta = get_theta(N, time, phase, R_mu, R_sigma, R_lambda, R_beta, x_mu, x_sigma, x_lambda, x_beta, x_alpha, delta_mu, delta_sigma, delta_lambda, delta_beta, delta_alpha);
}
model {

  // informative priors for intercepts, to improve identifiability
  delta_mu[1, 1] ~ normal(-2, .25);
  delta_sigma[1, 1] ~ normal(-2, .25);
  delta_lambda[1, 1] ~ normal(0, .25); 
  delta_beta[1, 1] ~ normal(1, .25);
  delta_alpha[1] ~ normal(0, .25);
  
  // priors for additional process's parameters, if applicable
  if (R_mu > 1) {
    delta_mu[1, R_mu] ~ normal(-1, .25);
  }
  if (R_sigma > 1) {
    delta_sigma[1, R_sigma] ~ normal(-1, .25);
  }
  if (R_lambda > 1) {
    delta_lambda[1, R_lambda] ~ normal(0, .25);
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
  if (K_lambda > 1) {
    to_vector(delta_lambda[2:K_lambda, ]) ~ normal(0, .5);
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
