data {
  int<lower=0> n;       // Number of years
  int<lower=0> C[n];    // Count
  vector[n] year;       // Year
}

parameters {
  real<lower=-20,upper=20> alpha;
  real<lower=-10,upper=10> beta1;
  real<lower=-10,upper=10> beta2;
  real<lower=-10,upper=10> beta3;
}

transformed parameters {
  vector[n] log_lambda;

  log_lambda <- alpha +
                beta1 * year +
                beta2 * year .* year +
                beta3 * year .* year .* year;
}

model {
  // Priors
  alpha ~ uniform(-20, 20);
  beta1 ~ uniform(-10, 10);
  beta2 ~ uniform(-10, 10);
  beta3 ~ uniform(-10, 10);

  // Likelihood
  C ~ poisson_log(log_lambda);
}

generated quantities {
  vector[n] lambda;

  lambda <- exp(log_lambda);
}
