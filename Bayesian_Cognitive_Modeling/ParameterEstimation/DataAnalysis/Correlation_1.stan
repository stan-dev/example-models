// #### Notes to Stan model #######################################################
// ## 1) Multivariate normal distribution in Stan uses covariance matrix instead of 
// ##    precision matrix.
// ## 2) Multivariate normal distribution can be (and is) also vectorized.
// ## 3) Warnings may occur during sampling, ignore them.
// ################################################################################

// Pearson Correlation
data { 
  int<lower=0> n;
  array[n] vector[2] x;
}

parameters {
  vector[2] mu;
  vector<lower=0>[2] lambda;
  real<lower=-1, upper=1> r;
}

transformed parameters {
  vector<lower=0>[2] sigma;
  cov_matrix[2] T;

  // Reparameterization
  sigma[1] = inv_sqrt(lambda[1]);
  sigma[2] = inv_sqrt(lambda[2]);
  
  T[1,1] = square(sigma[1]);
  T[1,2] = r * sigma[1] * sigma[2];
  T[2,1] = r * sigma[1] * sigma[2];
  T[2,2] = square(sigma[2]);
}

model {
  // Priors
  mu ~ normal(0, sqrt(1000));
  lambda ~ gamma(0.001, 0.001);
  
  // Data
  x ~ multi_normal(mu, T);
}
