// response preparation model with flexible inputs
// inputs determine number of processes and covariates for each param
functions {
  vector get_theta(int N, int[] id, vector time, vector x, vector mu0, vector sigma0, vector beta0, vector alpha0, real delta_mu, real delta_sigma) {
    
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
      alpha[n] = inv_logit(alpha0[id[n]]);

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
  real<lower=0> mu_scale;
  real<lower=0> sigma_scale;
  real<lower=0> beta_scale;
  real<lower=0> alpha_scale;
  real mu_loc;
  real sigma_loc;
  real beta_loc;
  real alpha_loc;
  real delta_mu;
  real delta_sigma;
  vector[J] mu_raw;
  vector[J] sigma_raw;
  vector[J] beta_raw;
  vector[J] alpha_raw;
}
transformed parameters {
  vector[N] theta;
  vector[J] mu0;
  vector[J] sigma0;
  vector[J] beta0;
  vector[J] alpha0;

  mu0 = mu_raw * mu_scale + mu_loc;
  sigma0 = sigma_raw * sigma_scale + sigma_loc;
  beta0 = beta_raw * beta_scale + beta_loc;
  alpha0 = alpha_raw * alpha_scale + alpha_loc;

  theta = get_theta(N, id, time, x, mu0, sigma0, beta0, alpha0, delta_mu, delta_sigma);
}
model {

  // priors on group-level locations and scales of intercepts
  mu_loc ~ normal(-.5, .5);
  sigma_loc ~ normal(-2, .5);
  beta_loc ~ normal(2, .5);
  alpha_loc ~ normal(-1, .5);

  mu_scale ~ normal(.5, .5);
  sigma_scale ~ normal(.5, .5);
  beta_scale ~ normal(.5, .5);
  alpha_scale ~ normal(.5, .5);

  // priors on *raw* id-level intercepts
  mu_raw ~ normal(0, 1); 
  sigma_raw ~ normal(0, 1);
  beta_raw ~ normal(0, 1);
  alpha_raw ~ normal(0, 1);

  // fixed effects 
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
