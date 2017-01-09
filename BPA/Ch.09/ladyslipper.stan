// -------------------------------------------------
// States (S):
// 1 vegetative
// 2 flowering
// 3 dormant
// 4 dead
// Observations (O):
// 1 seen vegetative
// 2 seen flowering
// 3 not seen
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
  simplex[3] po[4, nind, n_occ_minus_1];

  for (i in 1:nind) {
    first[i] = first_capture(y[i]);

    for (t in 1:n_occ_minus_1) {
      // Define probabilities of O(t) given S(t)
      po[1, i, t, 1] = 1.0;
      po[1, i, t, 2] = 0.0;
      po[1, i, t, 3] = 0.0;
      po[2, i, t, 1] = 0.0;
      po[2, i, t, 2] = 1.0;
      po[2, i, t, 3] = 0.0;
      po[3, i, t, 1] = 0.0;
      po[3, i, t, 2] = 0.0;
      po[3, i, t, 3] = 1.0;
      po[4, i, t, 1] = 0.0;
      po[4, i, t, 2] = 0.0;
      po[4, i, t, 3] = 1.0;
    }
  }
}

parameters {
  vector<lower=0,upper=1>[n_occasions-1] s; // Survival probability
  vector<lower=0>[3] a;
  vector<lower=0>[3] b;
  vector<lower=0>[3] c;
}

transformed parameters {
  simplex[3] psiD; // Transitions from dormant
  simplex[3] psiV; // Transitions from vegetative
  simplex[3] psiF; // Transitions from flowering
  simplex[4] ps[4, nind, n_occasions-1];

  // Constraints
  for (i in 1:3){
    psiD[i] = a[i] / sum(a);
    psiV[i] = b[i] / sum(b);
    psiF[i] = c[i] / sum(c);
  }

  // Define state-transition and observation matrices
  for (i in 1:nind) {
    // Define probabilities of state S(t+1) given S(t)
    for (t in 1:n_occ_minus_1) {
      ps[1, i, t, 1] = s[t] * psiV[1];
      ps[1, i, t, 2] = s[t] * psiV[2];
      ps[1, i, t, 3] = s[t] * psiV[3];
      ps[1, i, t, 4] = 1.0 - s[t];
      ps[2, i, t, 1] = s[t] * psiF[1];
      ps[2, i, t, 2] = s[t] * psiF[2];
      ps[2, i, t, 3] = s[t] * psiF[3];
      ps[2, i, t, 4] = 1.0 - s[t];
      ps[3, i, t, 1] = s[t] * psiD[1];
      ps[3, i, t, 2] = s[t] * psiD[2];
      ps[3, i, t, 3] = s[t] * psiD[3];
      ps[3, i, t, 4] = 1.0 - s[t];
      ps[4, i, t, 1] = 0.0;
      ps[4, i, t, 2] = 0.0;
      ps[4, i, t, 3] = 0.0;
      ps[4, i, t, 4] = 1.0;
    }
  }
}

model {
  real acc[4];
  vector[4] gamma[n_occasions];

  // Priors
  // Survival: uniform
  //  s ~ uniform(0, 1);

  // Transitions: gamma priors
  a ~ gamma(1, 1);
  b ~ gamma(1, 1);
  c ~ gamma(1, 1);

  // Likelihood
  // Forward algorithm derived from Stan Modeling Language
  // User's Guide and Reference Manual
  for (i in 1:nind) {
    if (first[i] > 0) {
      for (k in 1:4)
        gamma[first[i], k] = (y[i, first[i]] == k);

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
