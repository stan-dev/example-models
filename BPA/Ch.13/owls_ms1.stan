// Multistate occupancy models

data {
  int<lower=1> R;               // Number of sites
  int<lower=1> T;               // Number of temporal replications
  int<lower=0,upper=3> y[R, T]; // Observed data (0:NA)
}

parameters {
  real<lower=0,upper=1> p2;     // Detection prob. at a site w/o repro.
  real<lower=0,upper=1> psi;    // Occupancy prob.
  real<lower=0,upper=1> r;      // Reproduction prob.
  vector<lower=0>[3] beta;
}

transformed parameters {
  simplex[3] p3;                // Detectin prob.
  simplex[3] phi[R];            // State vector
  simplex[3] p[3, T];           // Observation matrix

  p3 = beta / sum(beta);       // Induce Dirichlet prior

  // Define state vector
  for (s in 1:R) {
    phi[s, 1] = 1 - psi;          // Prob. of non-occupation
    phi[s, 2] = psi * (1 - r);    // Prob. of occupancy without repro
    phi[s, 3] = psi * r;            // Prob. of occupancy and repro
  }

  // Define observation matrix
  // Order of indices: true state, time, observed state
  for (t in 1:T) {
    p[1, t, 1] = 1;
    p[1, t, 2] = 0;
    p[1, t, 3] = 0;
    p[2, t, 1] = 1 - p2;
    p[2, t, 2] = p2;
    p[2, t, 3] = 0;
    p[3, t, 1] = p3[1];
    p[3, t, 2] = p3[2];
    p[3, t, 3] = p3[3];
  }
}

model {
  // Priors
  // Flat priros are implicitly used on psi, r and p2.
  beta ~ gamma(1, 1); // Induce Dirichlet prior

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
