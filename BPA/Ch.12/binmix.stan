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
    max_y[i] <- max(y[i]);
}

parameters {
  real<lower=0> lambda;    // Mean population size
  real<lower=0,upper=1> p; // Detection probability
}

model {
  // Priors
  lambda ~ gamma(0.005, 0.005);
  // An implicit flat prior [0, 1] is used for p.

  // Likelihood
  for (i in 1:R) {
    vector[K+1] lp;

    for (n in max_y[i]:K) {
      lp[n + 1] <- poisson_log(n, lambda)
        + binomial_log(y[i], n, p);
    }
    increment_log_prob(log_sum_exp(lp[(max_y[i] + 1):(K + 1)]));
  }
}
