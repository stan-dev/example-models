data {
  int<lower=0> n;       // Number of years
  int<lower=0> C[n];    // Count
  vector[n] year;       // Year
}

transformed data {
  vector[n] year_squared;
  vector[n] year_cubed;

  year_squared = year .* year;
  year_cubed = year .* year .* year;
}

parameters {
  real<lower=-20,upper=20> alpha;
  real<lower=-10,upper=10> beta1;
  real<lower=-10,upper=10> beta2;
  real<lower=-10,upper=10> beta3;
}

transformed parameters {
  vector[n] log_lambda;

  log_lambda = alpha
             + beta1 * year +
             + beta2 * year_squared +
             + beta3 * year_cubed;
}

model {
  // Implicit uniform priors are used.

  // Likelihood
  C ~ poisson_log(log_lambda);
}

generated quantities {
  vector[n] lambda;

  lambda = exp(log_lambda);
}
