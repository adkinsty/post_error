// response preparation model with 2 processes and no covariates
functions {
    vector get_theta(int h, int g, int N, vector time, vector congr, vector mu, vector sigma, vector beta, real alpha) {

        vector[2] phi;
        vector[2] not_phi;
        vector[4] psi; 

        vector[N] theta; // p(corr | t)

        for (n in 1:N) {

            phi[h] = normal_cdf(time[n], mu[h], sigma[h]); //  p(h prep | t)
            phi[g] = normal_cdf(time[n], mu[g], sigma[g]); //  p(g prep | t)

            not_phi[h] = 1 - phi[h]; //  p(not h prep | t)
            not_phi[g] = 1 - phi[g]; //  p(not g prep | t)

            psi[1] = not_phi[h] * not_phi[g] * alpha;   // p(corr | not h prep & not g prep)
            
            if (congr[n] > 0) {
                psi[2] = phi[h]  * not_phi[g] * beta[h];  // p(corr | h prep & not g prep)
            } else {
                psi[2] = phi[h]  * not_phi[g] * (1 - beta[h]); 
            }
            
            psi[3] = not_phi[h] * phi[g] * beta[g]; // p(corr | not h prep & g prep)
            psi[4] = phi[h]  * phi[g] * beta[g];    // p(corr | h prep & g prep)

            theta[n] = sum(psi);  // response prob (@t)
        }
        return(theta);
    }
}

data {
  
  int<lower=0> N; // n observations
  int<lower=0,upper=1> y[N]; // observations
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
  ordered[2] mu_loc;
  vector[2] sigma_loc;
  vector[2] beta_loc;
  real alpha_loc;
}

transformed parameters {
  vector[N] theta;
  vector[2] mu;
  vector[2] sigma;
  vector[2] beta;
  real alpha;

  mu = inv_logit(mu_loc);
  sigma = inv_logit(sigma_loc);
  beta = inv_logit(beta_loc);  
  alpha = inv_logit(alpha_loc);
  
  // probability of correct reponse at each time given parameters
  theta = get_theta(h, g, N, time, congr, mu, sigma, beta, alpha);
}

model {

  // priors on group-level means of intercepts and slopes
  mu_loc ~ normal(-1, .5);
  sigma_loc ~ normal(-2, .5);
  beta_loc ~ normal(2, .5);
  alpha_loc ~ normal(0, .5); // assumes 2 buttons

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
