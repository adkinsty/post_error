// response preparation model with flexible inputs
functions {
  vector get_theta(int N, vector time, real mu_loc, real sigma_loc, real beta_loc, real alpha_loc) { //, vector delta_mu, vector delta_sigma, vector delta_beta) {
    
    vector[N] mu; // prep mean
    vector[N] sigma; // prep sd
    vector[N] beta; // p(corr | prep)
    vector[N] alpha; // p(corr | not prep)

    real phi; // p(prep | t)
    real not_phi; // p(not prep | t)
    vector[2] psi;  // p(corr | prep) & p(corr | not prep) 

    vector[N] theta; // p(corr | t)

    for (n in 1:N) {

      mu[n]    = inv_logit(mu_loc);
      sigma[n] = inv_logit(sigma_loc); 
      beta[n]  = inv_logit(beta_loc);
      alpha[n] = inv_logit(alpha_loc);

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
  int<lower=0,upper=1> y[N]; // observations
  vector<lower=0,upper=1>[N] time; // time in seconds

  int<lower=0,upper=1> gq; 
  int<lower=0,upper=1> prior_only;

}
parameters {

  real mu_loc;
  real sigma_loc;
  real beta_loc;
  real alpha_loc;

}
transformed parameters {
  vector[N] theta;

  theta = get_theta(N, time, mu_loc, sigma_loc, beta_loc, alpha_loc); //, delta_mu, delta_sigma, delta_beta);
}
model {

  // priors on group-level locations and scales of intercepts
  mu_loc ~ normal(-.5, .5);
  sigma_loc ~ normal(-2, .5);
  beta_loc ~ normal(2, .5);
  alpha_loc ~ normal(-1, .5);

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
