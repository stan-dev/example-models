// Jolly-Seber model as s multistate model

//--------------------------------------
// States (S):
// 1 not yet entered
// 2 alive
// 3 dead
// Observations (O):
// 1 seen
// 2 not seen
//--------------------------------------

functions {
  /**
   * Returns delta, where
   * delta[n] = gamma[n] * PROD_{m < n} (1 - gamma[m])
   * (Thanks to Dr. Carpenter)
   *
   * @param gamma Vector of probability sequence
   *
   * @return Vector of complementary probability sequence
   */
  vector seq_cprob(vector gamma) {
    int N = rows(gamma);
    vector[N] log_cprob;
    real log_residual_prob = 0;

    for (n in 1:N) {
      log_cprob[n] = log(gamma[n]) + log_residual_prob;
      log_residual_prob = log_residual_prob + log(1 - gamma[n]);
    }
    return exp(log_cprob);
  }
}

data {
  int<lower=0> M;                               // Augmented sample size
  int<lower=0> n_occasions;                     // Number of capture occasions
  int<lower=1,upper=2> y[M, n_occasions];       // Augmented capture-history
}

transformed data {
  int n_occ_minus_1 = n_occasions - 1;
}

parameters {
  vector<lower=0,upper=1>[n_occ_minus_1] gamma; // Removal entry probabilities
  real<lower=0,upper=1> mean_phi;               // Mean survival
  real<lower=0,upper=1> mean_p;                 // Mean capture
}

transformed parameters {
  vector<lower=0,upper=1>[n_occ_minus_1] phi; // Survival probability
  vector<lower=0,upper=1>[n_occ_minus_1] p;   // Capture probability
  simplex[3] ps[3, M, n_occ_minus_1];
  simplex[2] po[3, M, n_occ_minus_1];

  // Constraints
  for (t in 1:n_occ_minus_1) {
    phi[t] = mean_phi;
    p[t] = mean_p;
  }

  // Define state-transition and observation matrices
  for (i in 1:M) {
    // Define probabilities of state S(t+1) given S(t)
    for (t in 1:n_occ_minus_1) {
      ps[1, i, t, 1] = 1.0 - gamma[t];
      ps[1, i, t, 2] = gamma[t];
      ps[1, i, t, 3] = 0.0;
      ps[2, i, t, 1] = 0.0;
      ps[2, i, t, 2] = phi[t];
      ps[2, i, t, 3] = 1 - phi[t];
      ps[3, i, t, 1] = 0.0;
      ps[3, i, t, 2] = 0.0;
      ps[3, i, t, 3] = 1.0;

      // Define probabilities of O(t) given S(t)
      po[1, i, t, 1] = 0.0;
      po[1, i, t, 2] = 1.0;
      po[2, i, t, 1] = p[t];
      po[2, i, t, 2] = 1.0 - p[t];
      po[3, i, t, 1] = 0.0;
      po[3, i, t, 2] = 1.0;
    }
  }
}

model {
  real acc[3];
  vector[3] gam[n_occasions];

  // Priors
  // Uniform priors are implicitly defined.
  //  gamma ~ uniform(0, 1);
  //  mean_phi ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);

  // Likelihood
  // Forward algorithm derived from Stan Modeling Language
  // User's Guide and Reference Manual
  for (i in 1:M) {
    // Make sure that all M individuals are in state 1 at t=1
    gam[1, 1] = 1.0;
    gam[1, 2] = 0.0;
    gam[1, 3] = 0.0;

    for (t in 2:n_occasions) {
      for (k in 1:3) {
        for (j in 1:3)
          acc[j] = gam[t - 1, j] * ps[j, i, t - 1, k]
                  * po[k, i, t - 1, y[i, t]];
        gam[t, k] = sum(acc);
      }
    }
    target += log(sum(gam[n_occasions]));
  }
}

generated quantities {
  real<lower=0,upper=1> psi;                // Inclusion probability
  vector<lower=0,upper=1>[n_occ_minus_1] b; // Entry probability
  int<lower=0> Nsuper;                      // Superpopulation size
  int<lower=0> N[n_occ_minus_1];            // Actual population size
  int<lower=0> B[n_occ_minus_1];            // Number of entries
  int<lower=1,upper=3> z[M, n_occasions];   // Latent state

  // Generate z[]
  for (i in 1:M) {
    z[i, 1] = 1;
    for (t in 2:n_occasions)
      z[i, t] = categorical_rng(ps[z[i, t - 1], i, t - 1]);
  }

  // Calculate derived population parameters
  {
    vector[n_occ_minus_1] cprob  = seq_cprob(gamma[:n_occ_minus_1]);
    int al[M, n_occ_minus_1];
    int d[M, n_occ_minus_1];
    int alive[M];
    int w[M];

    psi = sum(cprob);
    b = cprob / psi;

    for (i in 1:M) {
      for (t in 2:n_occasions)
        al[i, t - 1] = (z[i, t] == 2);
      for (t in 1:n_occ_minus_1)
        d[i, t] = (z[i, t] == al[i, t]);
      alive[i] = sum(al[i]);
    }

    for (t in 1:n_occ_minus_1) {
      N[t] = sum(al[, t]);
      B[t] = sum(d[, t]);
    }
    for (i in 1:M)
      w[i] = 1 - !alive[i];
    Nsuper = sum(w);
  }
}
