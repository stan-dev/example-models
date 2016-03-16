// The simplest possible site-occupancy model

data {
  int<lower=1> R;                 // Number of sites
  int<lower=1> T;                 // Number of temporal replications
  int<lower=0,upper=1> y[R, T];   // Observation
}

transformed data {
  int<lower=0,upper=T> sum_y[R];  // Number of occupation for each site
  int<lower=0,upper=R> occ_obs;   // Number of observed occupied sites

  occ_obs <- 0;
  for (i in 1:R) {
    sum_y[i] <- sum(y[i]);
    if (sum_y[i])
      occ_obs <- occ_obs + 1;
  }
}

parameters {
  real<lower=0,upper=1> psi;      // Occupancy probability
  real<lower=0,upper=1> p;        // Detection probability
}

model {
  // Priors
  // Flat priors are implicitly used on psi and p.

  // Likelihood
  for (i in 1:R) {
    if (sum_y[i]) {    // Occurred and observed
      1 ~ bernoulli(psi);
      y[i] ~ bernoulli(p);
    } else {
                                     // Occurred and not observed
      increment_log_prob(log_sum_exp(bernoulli_log(1, psi)
                                     + bernoulli_log(0, p) * T,
                                     // Not occurred
                                     bernoulli_log(0, psi)));
    }
  }
}

generated quantities {
  real<lower=occ_obs> occ_fs;

  // Observed number of occupied sites
  // + Number of sites occupied and not observed
  occ_fs <- occ_obs + binomial_rng(R, psi * (1.0 - p)^T);
}
