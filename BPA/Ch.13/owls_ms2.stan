// Multistate occupancy models

data {
  int<lower=1> R;                // Number of sites
  int<lower=1> T;                // Number of temporal replications
  int<lower=0,upper=3> y[R, T];  // Observed data (0:NA)
}

parameters {
  vector<lower=0,upper=1>[T] p2; // Detection prob. at a site w/o repro.
  real<lower=0,upper=1> psi;     // Occupancy prob.
  real<lower=0,upper=1> r;       // Reproduction prob.
  vector<lower=0>[T] beta[3];
}

transformed parameters {
  simplex[3] p3[T];              // Detectin prob.
  simplex[3] phi[R];             // State vector
  simplex[3] p[3, T];            // Observation matrix

  for (t in 1:T)
    for (i in 1:3)
      p3[t, i] = beta[i, t] / sum(beta[, t]);

  // Define state vector
  phi[, 1] = rep_array(1 - psi, R);       // Prob. of non-occupation
  phi[, 2] = rep_array(psi * (1 - r), R); // Prob. of occupancy without repro
  phi[, 3] = rep_array(psi * r, R);       // Prob. of occupancy and repro

  // Define observation matrix
  // Order of indices: true state, time, observed state
  p[1, , 1] = rep_array(1, T);
  p[1, , 2] = rep_array(0, T);
  p[1, , 3] = rep_array(0, T);
  for (t in 1:T) {
    p[2, t, 1] = 1 - p2[t];
    p[2, t, 2] = p2[t];
  }
  p[2, , 3] = rep_array(0, T);
  p[3, , 1] = p3[, 1];
  p[3, , 2] = p3[, 2];
  p[3, , 3] = p3[, 3];
}

model {
  real acc[3];
  vector[3] gamma[T];

  // Priors
  // Flat priros are implicitly used on psi, r and p2.
  for (t in 1:T)
    beta[1:3, t] ~ gamma(1, 1); // Induce Dirichlet prior

  // Likelihood
  for (s in 1:R) {
    vector[3] lp;

    for (k in 1:3) {
      lp[k] = categorical_lpmf(k | phi[s]);
      for (t in 1:T) {
        if (y[s, t])
          lp[k] = lp[k] + categorical_lpmf(y[s, t] | p[k, t]);
      }
    }
    target += log_sum_exp(lp);
  }
}

generated quantities {
  int<lower=1,upper=3> z[R];    // State
  int occ[3, R];
  real n_occ[3];                // Number of each state

  for (s in 1:R)
    z[s] = categorical_rng(phi[s]);
  for (i in 1:3)
    for (s in 1:R)
      occ[i, s] = (z[s] == i);
  for (i in 1:3)
    n_occ[i] = sum(occ[i]);
}
