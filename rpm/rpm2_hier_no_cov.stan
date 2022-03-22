// response preparation model with 2 processes and no covariates
functions {
  vector get_theta(int h, int g, int N, int[] id, vector time, vector congr, matrix mu, matrix sigma, matrix beta, vector alpha) {
    
    vector[2] phi;
    vector[2] not_phi;
    vector[4] psi; 

    vector[N] theta; // p(corr | t)

    for (n in 1:N) {

      phi[h] = normal_cdf(time[n], mu[id[n], h], sigma[id[n], h]); //  p(h prep | t)
      phi[g] = normal_cdf(time[n], mu[id[n], g], sigma[id[n], g]); //  p(g prep | t)

      not_phi[h] = 1 - phi[h]; //  p(not h prep | t)
      not_phi[g] = 1 - phi[g]; //  p(not g prep | t)

      psi[1] = not_phi[h] * not_phi[g] * alpha[id[n]];   // p(corr | not h prep & not g prep)
      
      if (congr[n] > 0) {
        psi[2] = phi[h]  * not_phi[g] * beta[id[n], h];  // p(corr | h prep & not g prep)
      } else {
        psi[2] = phi[h]  * not_phi[g] * (1 - beta[id[n], h]); 
      }
      
      psi[3] = not_phi[h] * phi[g] * beta[id[n], g]; // p(corr | not h prep & g prep)
      psi[4] = phi[h]  * phi[g] * beta[id[n], g];    // p(corr | h prep & g prep)

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

  vector[N] congr; // habit and goal equal?

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
  real<lower=0> alpha_scale;

  ordered[2] mu_loc;
  vector[2] sigma_loc;
  vector[2] beta_loc;
  real alpha_loc;

  matrix[J, 2] mu_raw;
  matrix[J, 2] sigma_raw;
  matrix[J, 2] beta_raw;
  vector[J] alpha_raw;

}

transformed parameters {
  vector[N] theta;
  matrix[J, 2] mu;
  matrix[J, 2] sigma;
  matrix[J, 2] beta;
  vector[J] alpha;

  mu[, h] = inv_logit(mu_raw[, h] * mu_scale[h] + mu_loc[h]);
  mu[, g] = inv_logit(mu_raw[, g] * mu_scale[g] + mu_loc[g]);

  sigma[, h] = inv_logit(sigma_raw[, h] * sigma_scale[h] + sigma_loc[h]);
  sigma[, g] = inv_logit(sigma_raw[, g] * sigma_scale[g] + sigma_loc[g]);  
  
  beta[, h] = inv_logit(beta_raw[, h] * beta_scale[h] + beta_loc[h]);
  beta[, g] = inv_logit(beta_raw[, g] * beta_scale[g] + beta_loc[g]);    
  
  alpha = inv_logit(alpha_raw * alpha_scale + alpha_loc);
  
  // probability of correct reponse at each time given parameters
  theta = get_theta(h, g, N, id, time, congr, mu, sigma, beta, alpha);
}

model {

  // priors on group-level sds of intercepts and slopes
  mu_scale ~ normal(0, .5);
  sigma_scale ~ normal(0, .5);
  beta_scale ~ normal(0, .5);
  alpha_scale ~ normal(0, .5);

  // priors on group-level means of intercepts and slopes
  mu_loc ~ normal(-1, .5);
  sigma_loc ~ normal(-2, .5);
  beta_loc ~ normal(2, .5);
  alpha_loc ~ normal(0, .5); // assumes 2 buttons

  // priors on *raw* id-level intercepts and slopes
  to_vector(mu_raw) ~ normal(0, 1); 
  to_vector(sigma_raw) ~ normal(0, 1);
  to_vector(beta_raw) ~ normal(0, 1);
  to_vector(alpha_raw) ~ normal(0, 1);

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
