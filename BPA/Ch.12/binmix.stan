// Simplest binomial mixture model
data {
  int<lower=0> R;       // Number of sites
  int<lower=0> T;       // Number of temporal replications
  int<lower=0> y[R, T]; // Counts
  int<lower=0> K;       // Upper bound of population size
}

transformed data {
  int<lower=0> max_y[R];

  for (i in 1:R)
    max_y[i] = max(y[i]);
}

parameters {
  real<lower=0> lambda;    // Mean population size
  real<lower=0,upper=1> p; // Detection probability
}

model {
  // Priors
  // A half Cauchy prior is used on lambda, instead of
  // gamma(0.005, 0.005) or uniform(0, 10) proposed in the book.
  lambda ~ cauchy(0, 10);
  // A flat prior [0, 1] is implicitly used on p.

  // Likelihood
  for (i in 1:R) {
    vector[K - max_y[i] + 1] lp;

    for (j in 1:(K - max_y[i] + 1))
      lp[j] = poisson_lpmf(max_y[i] + j - 1 | lambda)
             + binomial_lpmf(y[i] | max_y[i] + j - 1, p);
    target += log_sum_exp(lp);
  }
}
