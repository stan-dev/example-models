// Single-season occupancy model
data {
  int<lower=1> R;
  int<lower=1> T;
  int<lower=0,upper=1> y[R, T];
  int<lower=0,upper=1> edge[R];
  matrix[R, T] DATES;
  matrix[R, T] HOURS;
  int last[R];
}

transformed data {
  int<lower=0,upper=T> sum_y[R];
  int<lower=0,upper=R> occ_obs;  // Number of observed occupied sites
  matrix[R, T] DATES2;
  matrix[R, T] HOURS2;

  occ_obs = 0;
  for (i in 1:R) {
    sum_y[i] = sum(y[i]);
    if (sum_y[i])
      occ_obs = occ_obs + 1;
  }
  DATES2 = DATES .* DATES;
  HOURS2 = HOURS .* HOURS;
}

parameters {
  real alpha_psi;
  real beta_psi;
  real alpha_p;
  real beta1_p;
  real beta2_p;
  real beta3_p;
  real beta4_p;
}

transformed parameters {
  vector[R] logit_psi;  // Logit occupancy prob.
  matrix[R, T] logit_p; // Logit detection prob.

  for (i in 1:R)
    logit_psi[i] = alpha_psi + beta_psi * edge[i];
  logit_p = alpha_p
      + beta1_p * DATES + beta2_p * DATES2
      + beta3_p * HOURS + beta4_p * HOURS2;
}

model {
  // Priors
  alpha_psi ~ normal(0, 10);
  beta_psi ~ normal(0, 10);
  alpha_p ~ normal(0, 10);
  beta1_p ~ normal(0, 10);
  beta2_p ~ normal(0, 10);
  beta3_p ~ normal(0, 10);
  beta4_p ~ normal(0, 10);

  // Likelihood
  for (i in 1:R) {
    if (sum_y[i]) { // Occupied and observed
      target += bernoulli_logit_lpmf(1 |  logit_psi[i])
        + bernoulli_logit_lpmf(y[i, 1:last[i]] | logit_p[i, 1:last[i]]);
    } else {        // Never observed
                            // Occupied and not observed
      target += log_sum_exp(bernoulli_logit_lpmf(1 | logit_psi[i])
                            + bernoulli_logit_lpmf(0 | logit_p[i, 1:last[i]]),
                            // Not occupied
                            bernoulli_logit_lpmf(0 | logit_psi[i]));
    }
  }
}

generated quantities {
  real<lower=0,upper=1> mean_p = inv_logit(alpha_p);
  int occ_fs;       // Number of occupied sites
  real psi_con[R];  // prob present | data
  int z[R];         // occupancy indicator, 0/1
  
  for (i in 1:R) {
    if (sum_y[i] == 0) {  // species not detected
      real psi = inv_logit(logit_psi[i]);
      vector[last[i]] q = inv_logit(-logit_p[i, 1:last[i]])';  // q = 1 - p
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
