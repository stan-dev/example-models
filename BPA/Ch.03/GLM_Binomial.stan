data {
  int<lower=0> nyears; // Number of Years
  array[nyears] int<lower=0> C; // Counts
  array[nyears] int<lower=0> N; // Binomial Totals
  vector[nyears] year; // Year covariates
}
transformed data {
  vector[nyears] year_squared;
  
  year_squared = year .* year;
}
parameters {
  real alpha;
  real beta1;
  real beta2;
}
transformed parameters {
  vector[nyears] logit_p;
  
  // Linear predictor
  logit_p = alpha + beta1 * year + beta2 * year_squared;
}
model {
  // Priors
  alpha ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  beta2 ~ normal(0, 100);
  
  // Likelihood
  // Distribution for random part
  C ~ binomial_logit(N, logit_p);
}
generated quantities {
  array[nyears] real<lower=0, upper=1> p;
  
  for (i in 1 : nyears) {
    p[i] = inv_logit(logit_p[i]);
  }
}
