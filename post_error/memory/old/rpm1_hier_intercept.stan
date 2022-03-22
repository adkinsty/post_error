// response preparation model with flexible inputs
// inputs determine number of processes and covariates for each param
functions {
  vector get_theta(int N, int[] id, vector time, vector x, vector mu0, vector sigma0, vector beta0, real alpha0, real delta_mu, real delta_sigma) {
    
    vector[N] mu; // prep mean
    vector[N] sigma; // prep sd
    vector[N] beta; // p(corr | prep)
    vector[N] alpha; // p(corr | not prep)

    real phi; // p(prep | t)
    real not_phi; // p(not prep | t)
    vector[2] psi;  // p(corr | prep) & p(corr | not prep) 

    vector[N] theta; // p(corr | t)

    for (n in 1:N) {

      mu[n]    = inv_logit(mu0[id[n]] + delta_mu * x[n]);
      sigma[n] = inv_logit(sigma0[id[n]] + delta_sigma * x[n]);
      beta[n]  = inv_logit(beta0[id[n]]);
      alpha[n] = inv_logit(alpha0);

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
  
  int<lower=0> N; // n observations
  int<lower=0> J; // n subjects
  int<lower=0,upper=1> y[N]; // observations
  int<lower=0> id[N]; // subject ids
  vector<lower=0,upper=1>[N] time; // time in seconds

  vector[N] x; //covariate

  int<lower=0,upper=1> gq; 
  int<lower=0,upper=1> prior_only;

}
parameters {
  real<lower=0> mu_sd;
  real<lower=0> sigma_sd;
  real<lower=0> beta_sd;
  real mu_mean;
  real sigma_mean;
  real beta_mean;
  real alpha0;
  real delta_mu;
  real delta_sigma;
  vector[J] mu0;
  vector[J] sigma0;
  vector[J] beta0;
}
transformed parameters {
  vector[N] theta;
  theta = get_theta(N, id, time, x, mu0, sigma0, beta0, alpha0, delta_mu, delta_sigma);
}
model {

  // priors on group-level means and sds of intercepts
  mu_mean ~ normal(-.5, .5);
  sigma_mean ~ normal(-2, .5);
  beta_mean ~ normal(2, .5);
  mu_sd ~ normal(.5, .5);
  sigma_sd ~ normal(.5, .5);
  beta_sd ~ normal(.5, .5);

  // priors on id-level intercepts
  mu0 ~ normal(mu_mean, mu_sd); 
  sigma0 ~ normal(sigma_mean, sigma_sd);
  beta0 ~ normal(beta_mean, beta_sd);

  // fixed effects 
  alpha0 ~ normal(-1, .5);
  delta_mu ~ normal(0, .5);
  delta_sigma ~ normal(0, .5);

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
