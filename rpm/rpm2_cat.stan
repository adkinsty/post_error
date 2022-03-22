// response preparation model with 2 processes, hiearchical priors, and categorical responses
functions {
  vector[] get_theta(int h, int g, int N, vector time,  vector mu_loc, vector sigma_loc, vector beta_loc, real alpha_loc) { 
    
    matrix[N, 2] mu; // prep mean
    matrix[N, 2] sigma; // prep sd
    matrix[N, 2] beta; // p(corr | prep)
    vector[N] alpha; // p(corr | not prep)

    vector[2] cdf;
    matrix[4, 1] phi;
    matrix[3, 4] psi; 

    vector[3] theta[N]; // p(resp | t)

    for (n in 1:N) {

      mu[n, h]    = inv_logit(mu_loc[h]); 
      mu[n, g]    = inv_logit(mu_loc[g]); 

      sigma[n, h] = inv_logit(sigma_loc[h]); 
      sigma[n, g] = inv_logit(sigma_loc[g]); 

      beta[n, h] = inv_logit(beta_loc[h]); 
      beta[n, g] = inv_logit(beta_loc[g]);

      alpha[n] = inv_logit(alpha_loc);

      cdf[h] = normal_cdf(time[n], mu[n, h], sigma[n, h]); //  p(h prep | t)
      cdf[g] = normal_cdf(time[n], mu[n, g], sigma[n, g]); //  p(g prep | t)

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

      theta[n] = to_vector(psi * phi);


    }
    return(theta);
  }
}

data {
  
  int<lower=0> N; // n observations
  int<lower=1,upper=3> y[N]; // observations
  vector<lower=0,upper=1>[N] time; // time in seconds

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
  ordered[2] mu_loc;
  vector[2] sigma_loc;
  vector[2] beta_loc;
  real alpha_loc;
}
transformed parameters {

  vector[3] theta[N];
  vector[N] theta1;
  vector[N] theta2;
  vector[N] theta3;
  
  theta = get_theta(h, g, N, time, mu_loc, sigma_loc, beta_loc, alpha_loc);
  theta1 = to_vector(theta[ ,1]);
  theta2 = to_vector(theta[ ,2]);
  theta3 = to_vector(theta[ ,3]);
}

model {

  // priors on group-level intercepts
  mu_loc ~ normal(-1, .5);
  sigma_loc ~ normal(-2, .5);
  beta_loc ~ normal(2, .5);
  alpha_loc ~ normal(-2, .5);

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
