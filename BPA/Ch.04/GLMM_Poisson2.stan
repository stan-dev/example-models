data {
  int<lower=0> nsite; // Number of populations
  int<lower=0> nyear; // Number of years
  array[nyear, nsite] int<lower=0> C; // Counts
  vector[nyear] year; // Year covariate
}
transformed data {
  vector[nyear] year_squared;
  vector[nyear] year_cubed;
  
  year_squared = year .* year;
  year_cubed = year .* year .* year;
}
parameters {
  real mu;
  vector[nsite] alpha;
  array[nyear] real eps;
  array[3] real beta;
  real<lower=0, upper=2> sd_alpha;
  real<lower=0, upper=1> sd_year;
}
transformed parameters {
  array[nyear] vector[nsite] log_lambda;
  
  // Linear predictor including random site and random year effects
  for (i in 1 : nyear) {
    log_lambda[i] = alpha + beta[1] * year[i] + beta[2] * year_squared[i]
                    + beta[3] * year_cubed[i] + eps[i];
  }
}
model {
  // Priors
  
  // Random site effects
  alpha ~ normal(mu, sd_alpha);
  
  // Hyperparameter 1
  mu ~ normal(0, 10);
  
  // Hyperparameter 2
  //  sd_alpha ~ uniform(0, 2); // Implicitly defined
  
  beta ~ normal(0, 10);
  
  // Hyperparameter 3
  //  sd_year ~ uniform(0, 1); // Implicitly defined
  
  // Random year effects
  eps ~ normal(0, sd_year);
  
  // Likelihood
  for (i in 1 : nyear) {
    // Distribution for random part
    // Link function
    C[i] ~ poisson_log(log_lambda[i]);
  }
}
