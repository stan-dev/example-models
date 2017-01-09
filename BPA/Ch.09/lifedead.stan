// -------------------------------------------------
// States (S):
// 1 alive in study area
// 2 alive outside study area
// 3 recently dead and recovered
// 4 recently dead, but not recovered, or dead (absorbing)
// Observations (O):
// 1 seen alive
// 2 recovered dead
// 3 neither seen nor recovered
// -------------------------------------------------

functions {
  /**
   * Return an integer value denoting occasion of first capture.
   * This function is derived from Stan Modeling Language
   * User's Guide and Reference Manual.
   *
   * @param y         Observed values
   * @return Occasion of first capture
   */
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k] != 3)
        return k;
    return 0;
  }
}

data {
  int<lower=0> nind;
  int<lower=0> n_occasions;
  int<lower=1,upper=3> y[nind, n_occasions];
}

transformed data {
  int n_occ_minus_1 = n_occasions - 1;
  int<lower=0,upper=n_occasions> first[nind];

  for (i in 1:nind)
    first[i] = first_capture(y[i]);
}

parameters {
  real<lower=0,upper=1> mean_s;    // Mean survival
  real<lower=0,upper=1> mean_f;    // Mean fidelity
  real<lower=0,upper=1> mean_r;    // Mean recovery
  real<lower=0,upper=1> mean_p;    // Mean recapture
}

transformed parameters {
  vector<lower=0,upper=1>[n_occ_minus_1] s; // True survival probability
  vector<lower=0,upper=1>[n_occ_minus_1] F; // Fidelity probability
  vector<lower=0,upper=1>[n_occ_minus_1] r; // Recovery probability
  vector<lower=0,upper=1>[n_occ_minus_1] p; // Recapture/resighting probability
  simplex[4] ps[4, nind, n_occ_minus_1];
  simplex[3] po[4, nind, n_occ_minus_1];

  // Constraints
  for (t in 1:n_occ_minus_1) {
    s[t] = mean_s;
    F[t] = mean_f;
    r[t] = mean_r;
    p[t] = mean_p;
  }

  // Define state-transition and observation matrices
  for (i in 1:nind) {
    // Define probabilities of state S(t+1) given S(t)
    for (t in 1:n_occ_minus_1) {
      ps[1, i, t, 1] = s[t] * F[t];
      ps[1, i, t, 2] = s[t] * (1.0 - F[t]);
      ps[1, i, t, 3] = (1.0 - s[t]) * r[t];
      ps[1, i, t, 4] = (1.0 - s[t]) * (1.0 - r[t]);
      ps[2, i, t, 1] = 0.0;
      ps[2, i, t, 2] = s[t];
      ps[2, i, t, 3] = (1.0 - s[t]) * r[t];
      ps[2, i, t, 4] = (1.0 - s[t]) * (1.0 - r[t]);
      ps[3, i, t, 1] = 0.0;
      ps[3, i, t, 2] = 0.0;
      ps[3, i, t, 3] = 0.0;
      ps[3, i, t, 4] = 1.0;
      ps[4, i, t, 1] = 0.0;
      ps[4, i, t, 2] = 0.0;
      ps[4, i, t, 3] = 0.0;
      ps[4, i, t, 4] = 1.0;

      // Define probabilities of O(t) given S(t)
      po[1, i, t, 1] = p[t];
      po[1, i, t, 2] = 0.0;
      po[1, i, t, 3] = 1.0 - p[t];
      po[2, i, t, 1] = 0.0;
      po[2, i, t, 2] = 0.0;
      po[2, i, t, 3] = 1.0;
      po[3, i, t, 1] = 0.0;
      po[3, i, t, 2] = 1.0;
      po[3, i, t, 3] = 0.0;
      po[4, i, t, 1] = 0.0;
      po[4, i, t, 2] = 0.0;
      po[4, i, t, 3] = 1.0;
      }
   }
}

model {
  real acc[4];
  vector[4] gamma[n_occasions];

  // Uniform priors are implicitly defined.
  //  mean_s ~ uniform(0, 1);
  //  mean_f ~ uniform(0, 1);
  //  mean_r ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);

  // Likelihood
  // Forward algorithm derived from Stan Modeling Language
  // User's Guide and Reference Manual
  for (i in 1:nind) {
    if (first[i] > 0) {
      for (k in 1:4)
        gamma[first[i], k] = (k == y[i, first[i]]);

      for (t in (first[i] + 1):n_occasions) {
        for (k in 1:4) {
          for (j in 1:4)
            acc[j] = gamma[t - 1, j] * ps[j, i, t - 1, k]
                    * po[k, i, t - 1, y[i, t]];
          gamma[t, k] = sum(acc);
        }
      }
      target += log(sum(gamma[n_occasions]));
    }
  }
}
