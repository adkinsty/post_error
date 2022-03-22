// response preparation model with 2 processes, hiearchical priors, and categorical responses
functions {
  vector[] get_theta(int h, int g, int N, int J, int[] id, vector time,  matrix mu0, matrix sigma0, matrix beta0, vector alpha0) { //,vector x, vector congr, matrix delta_mu, matrix delta_sigma, matrix delta_beta) {
    
    matrix[N, 2] mu; // prep mean
    matrix[N, 2] sigma; // prep sd
    matrix[N, 2] beta; // p(corr | prep)
    vector[N] alpha; // p(corr | not prep)

    vector[2] cdf;
    matrix[4, 1] phi;
    matrix[3, 4] psi; 

    vector[3] theta[N]; // p(resp | t)

    for (n in 1:N) {

      mu[n, h]    = inv_logit(mu0[id[n], h]); // + delta_mu[id[n], h] * x[n]);
      mu[n, g]    = inv_logit(mu0[id[n], g]); // + delta_mu[id[n], g] * x[n]);

      sigma[n, h] = inv_logit(sigma0[id[n], h]); // + delta_sigma[id[n], h] * x[n]);
      sigma[n, g] = inv_logit(sigma0[id[n], g]); // + delta_sigma[id[n], g] * x[n]);

      beta[n, h] = inv_logit(beta0[id[n], h]); // + delta_beta[id[n], h] * x[n]);
      beta[n, g] = inv_logit(beta0[id[n], g]); // + delta_beta[id[n], g] * x[n]);

      alpha[n] = inv_logit(alpha0[id[n]]);

      cdf[h] = normal_cdf(time[n], mu[n, h], sigma[n, h]); //  p(h prep | t)
      cdf[g] = normal_cdf(time[n], mu[n, g], sigma[n, g]); //  p(g prep | t)


      // if (congr[n] > 0) {
        
      //   theta[n, 1] = cdf[h] * beta[n, h];
      //   theta[n, 2] = 0
      //   theta[n, 3] = (1 - cdf[h]) * alpha[n];


      // } else {
      phi[1, 1] = (1 - cdf[h]) * (1 - cdf[g]);
      phi[2, 1] = cdf[h] * (1 - cdf[g]);
      phi[3, 1] = (1 - cdf[h]) * cdf[g];
      phi[4, 1] = cdf[h] * cdf[g];
    
      psi[1, 1] = alpha[n];
      psi[1, 2] = (1 - beta[n, h]) / 3;
      psi[1, 3] = beta[n, g];
      psi[1, 4] = beta[n, g];

      psi[2, 1] = alpha[n];
      psi[2, 2] = beta[n, h];
      psi[2, 3] = (1 - beta[n, g]) / 3;
      psi[2, 4] = (1 - beta[n, g]) / 3;

      psi[3, 1] = 1 - 2 * alpha[n];
      psi[3, 2] = (2 * (1 - beta[n, h])) / 3;
      psi[3, 3] = (2 * (1 - beta[n, g])) / 3;
      psi[3, 4] = (2 * (1 - beta[n, g])) / 3;

      // }

      theta[n] = to_vector(psi * phi);
      // print(phi);
      // print(psi);
      // print(theta[n]);

    }
    return(theta);
  }
}

data {
  
  int<lower=0> N; // n observations
  int<lower=0> J; // n subjects
  int<lower=1,upper=3> y[N]; // observations
  int<lower=0> id[N]; // subject ids
  vector<lower=0,upper=1>[N] time; // time in seconds

  // vector[N] x; //covariate
  // vector[N] congr; // habit and goal equal?

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
  vector<lower=0>[2] mu_scale;
  vector<lower=0>[2] sigma_scale;
  vector<lower=0>[2] beta_scale;
  // vector<lower=0>[2] delta_mu_scale;
  // vector<lower=0>[2] delta_sigma_scale;
  // vector<lower=0>[2] delta_beta_scale;
  real<lower=0> alpha_scale;

  ordered[2] mu_loc;
  ordered[2] sigma_loc;
  ordered[2] beta_loc;
  // vector[2] delta_mu_loc;
  // vector[2] delta_sigma_loc;
  // vector[2] delta_beta_loc;
  real alpha_loc;

  matrix[J, 2] mu_raw;
  matrix[J, 2] sigma_raw;
  matrix[J, 2] beta_raw;
  // matrix[J, 2] delta_mu_raw;
  // matrix[J, 2] delta_sigma_raw;
  // matrix[J, 2] delta_beta_raw;
  vector[J] alpha_raw;

}
transformed parameters {
  matrix[J, 2] mu0;
  matrix[J, 2] sigma0;
  matrix[J, 2] beta0;
  // matrix[J, 2] delta_mu;
  // matrix[J, 2] delta_sigma;
  // matrix[J, 2] delta_beta;
  vector[J] alpha0;

  vector[3] theta[N];

  // id level intercepts
  mu0[, h] = mu_raw[, h] * mu_scale[h] + mu_loc[h];
  mu0[, g] = mu_raw[, g] * mu_scale[g] + mu_loc[g];

  sigma0[, h] = sigma_raw[, h] * sigma_scale[h] + sigma_loc[h];
  sigma0[, g] = sigma_raw[, g] * sigma_scale[g] + sigma_loc[g];  
  
  beta0[, h] = beta_raw[, h] * beta_scale[h] + beta_loc[h];
  beta0[, g] = beta_raw[, g] * beta_scale[g] + beta_loc[g];    
  
  alpha0 = alpha_raw * alpha_scale + alpha_loc;

  // id level slopes
  // delta_mu[, h] = delta_mu_raw[, h] * delta_mu_scale[h] + delta_mu_loc[h]
  // delta_mu[, g] = delta_mu_raw[, g] * delta_mu_scale[g] + delta_mu_loc[g];

  // delta_sigma[, h] = delta_sigma_raw[, h] * delta_sigma_scale[h] + delta_sigma_loc[h];
  // delta_sigma[, g] = delta_sigma_raw[, g] * delta_sigma_scale[g] + delta_sigma_loc[g];

  // delta_beta[, h] = delta_beta_raw[, h] * delta_beta_scale[h] + delta_beta_loc[h];
  // delta_beta[, g] = delta_beta_raw[, g] * delta_beta_scale[g] + delta_beta_loc[g];
  
  // response probabilities at each time given parameters
  theta = get_theta(h, g, N, J, id, time, mu0, sigma0, beta0, alpha0); //  x, congr, delta_mu, delta_sigma, delta_beta);
}

model {

  // priors on group-level sds of intercepts and slopes
  mu_scale ~ normal(0, .5);
  sigma_scale ~ normal(0, .5);
  beta_scale ~ normal(0, .5);
  alpha_scale ~ normal(0, .5);
  // delta_mu_scale ~ normal(0, .5);
  // delta_sigma_scale ~ normal(0, .5);
  // delta_beta_scale ~ normal(0, .5);

  // priors on group-level means of intercepts and slopes
  mu_loc ~ normal(-1, .5);
  sigma_loc ~ normal(-2, .5);
  beta_loc ~ normal(2, .5);
  alpha_loc ~ normal(0, .5); // assumes 2 buttons
  // delta_mu_loc ~ normal(0, .5);
  // delta_sigma_loc ~ normal(0, .5);
  // delta_beta_loc ~ normal(0, .5);

  // priors on *raw* id-level intercepts and slopes
  to_vector(mu_raw) ~ normal(0, 1); 
  to_vector(sigma_raw) ~ normal(0, 1);
  to_vector(beta_raw) ~ normal(0, 1);
  to_vector(alpha_raw) ~ normal(0, 1);
  // to_vector(delta_mu_raw) ~ normal(0, 1);
  // to_vector(delta_sigma_raw) ~ normal(0, 1);
  // to_vector(delta_beta_raw) ~ normal(0, 1);

  if (prior_only == 0) {
    for (n in 1:N) {
      target += categorical_lpmf(y[n] | theta[n]);
    }
  }

}

generated quantities {

  vector[N * gq] log_lik;

  if (gq == 1) {

    for (n in 1:N) {
      
      log_lik[n] = categorical_lpmf(y[n] | theta[n]);

    }
  }
}
