// The simplest possible site-occupancy model

data {
  int<lower=1> R;                 // Number of sites
  int<lower=1> T;                 // Number of temporal replications
  int<lower=0,upper=1> y[R, T];   // Observation
}

transformed data {
  int<lower=0,upper=T> sum_y[R];  // Number of detections for each site
  int<lower=0,upper=R> occ_obs;   // Number of observed occupied sites

  occ_obs = 0;
  for (i in 1:R) {
    sum_y[i] = sum(y[i]);
    if (sum_y[i])
      occ_obs = occ_obs + 1;
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
      target += log_sum_exp(bernoulli_lpmf(1 | psi)
                            + bernoulli_lpmf(0 | p) * T,
                            // Not occurred
                            bernoulli_lpmf(0 | psi));
    }
  }
}

generated quantities {
  int<lower=occ_obs, upper=R> occ_fs;
  real psi_nd;  // prob occurred given not detected
  
  psi_nd = (psi * (1 - p)^T) / (psi * (1 - p)^T + (1 - psi));
  occ_fs = occ_obs + binomial_rng(R - occ_obs, psi_nd);
}
