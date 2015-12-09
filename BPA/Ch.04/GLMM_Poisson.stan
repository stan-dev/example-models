data {
  int<lower=0> n;       // Number of years
  int<lower=0> C[n];    // Counts
  real year[n];         // Year
}

parameters {
  real alpha;
  real beta1;
  real beta2;
  real beta3;
  real eps[n];          // Year effects
  real<lower=0> sigma;
}

transformed parameters {
  real log_lambda[n];

  // Linear predictor incl. random year effect
 for (i in 1:n)
   log_lambda[i] <- alpha +
                    beta1 * year[i] +
                    beta2 * pow(year[i], 2) +
                    beta3 * pow(year[i], 3) +
                    eps[i];
}

model {
  // Priors
  alpha ~ uniform(-20, 20);
  beta1 ~ uniform(-10, 10);
  beta2 ~ uniform(-10, 10);
  beta3 ~ uniform(-10, 10);
  sigma ~ uniform(0, 5);

  // Likelihood
  for (i in 1:n) {
    C[i] ~ poisson_log(log_lambda[i]);
    eps[i] ~ normal(0, sigma);
   }
}

generated quantities {
  real<lower=0> lambda[n];

  for (i in 1:n)
    lambda[i] <- exp(log_lambda[i]);
}
