// Site-occupancy models with covariates

data {
  int<lower=1> R;                 // Number of sites
  int<lower=1> T;                 // Number of temporal replications
  int<lower=0,upper=1> y[R, T];   // Observation
  vector[R] X;                    // Covariate
}

transformed data {
  int<lower=0,upper=T> sum_y[R];  // Number of occupation for each site
  int<lower=0,upper=R> occ_obs;   // Number of observed occupied sites

  occ_obs = 0;
  for (i in 1:R) {
    sum_y[i] = sum(y[i]);
    if (sum_y[i])
      occ_obs = occ_obs + 1;
  }
}

parameters {
  real alpha_occ;
  real beta_occ;
  real alpha_p;
  real beta_p;
}

transformed parameters {
  vector[R] logit_psi;            // Logit occupancy probability
  matrix[R, T] logit_p;           // Logit detection probability

  logit_psi = alpha_occ + beta_occ * X;
  logit_p = rep_matrix(alpha_p + beta_p * X, T);
}

model {
  // Priors
  // Improper flat priors are implicitly used on
  // alpha_occ, beta_occ, alpha_p and beta_p.

  // Likelihood
  for (i in 1:R) {
    if (sum_y[i]) {    // Occurred and observed
      1 ~ bernoulli_logit(logit_psi[i]);
      y[i] ~ bernoulli_logit(logit_p[i]);
    } else {
                            // Occurred and not observed
      target += log_sum_exp(bernoulli_logit_lpmf(1 | logit_psi[i])
                            + bernoulli_logit_lpmf(0 | logit_p[i]),
                            // Not occurred
                            bernoulli_logit_lpmf(0 | logit_psi[i]));
    }
  }
}

generated quantities {
  int occ_fs;       // Number of occupied sites
  real psi_con[R];  // prob occupied conditional on data
  int z[R];         // occupancy indicator, 0/1
  
  for (i in 1:R) {
    if (sum_y[i] == 0) {  // species not detected
      real psi = inv_logit(logit_psi[i]);
      vector[T] q = inv_logit(-logit_p[i])';  // q = 1 - p
      real qT = prod(q[]);
      psi_con[i] = (psi * qT) / (psi * qT + (1 - psi));
      z[i] = bernoulli_rng(psi_con[i]);
    } else {             // species detected at least once
      psi_con[i] = 1;
      z[i] = 1;
    }
  }
  occ_fs = sum(z);
}
