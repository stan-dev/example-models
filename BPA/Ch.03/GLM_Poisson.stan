data {
  int<lower=0> n;       // Number of years
  int<lower=0> C[n];    // Count
  real year[n];         // Year
}

parameters {
  real alpha;
  real beta1;
  real beta2;
  real beta3;
}

transformed parameters {
  vector[n] log_lambda;

  for (i in 1:n)
    log_lambda[i] <- alpha +
                     beta1 * year[i] +
                     beta2 * pow(year[i], 2) +
                     beta3 * pow(year[i], 3);
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
