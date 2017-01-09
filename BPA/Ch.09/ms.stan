// -------------------------------------------------
// States (S):
// 1 alive at A
// 2 alive at B
// 3 dead
// Observations (O):
// 1 seen at A
// 2 seen at B
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

  for (i in 1:nind)
    first[i] = first_capture(y[i]);
}

parameters {
  real<lower=0,upper=1> mean_phi[2]; // Mean state-spec. survival
  real<lower=0,upper=1> mean_psi[2]; // Mean transitions
  real<lower=0,upper=1> mean_p[2];   // Mean state-spec. recapture
}

transformed parameters {
  vector<lower=0,upper=1>[n_occ_minus_1] phiA;  // Survival prob. at site A
  vector<lower=0,upper=1>[n_occ_minus_1] phiB;  // Survival prob. at site B
  vector<lower=0,upper=1>[n_occ_minus_1] psiAB; // Movement prob. from site A to site B
  vector<lower=0,upper=1>[n_occ_minus_1] psiBA; // Movement prob. from site B to site A
  vector<lower=0,upper=1>[n_occ_minus_1] pA;    // Recapture prob. at site A
  vector<lower=0,upper=1>[n_occ_minus_1] pB;    // Recapture prob. at site B
  simplex[3] ps[3, nind, n_occ_minus_1];
  simplex[3] po[3, nind, n_occ_minus_1];

  // Constraints
  for (t in 1:n_occ_minus_1) {
    phiA[t] = mean_phi[1];
    phiB[t] = mean_phi[2];
    psiAB[t] = mean_psi[1];
    psiBA[t] = mean_psi[2];
    pA[t] = mean_p[1];
    pB[t] = mean_p[2];
  }

  // Define state-transition and observation matrices
  for (i in 1:nind) {
    // Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n_occasions - 1)) {
      ps[1, i, t, 1] = phiA[t] * (1.0 - psiAB[t]);
      ps[1, i, t, 2] = phiA[t] * psiAB[t];
      ps[1, i, t, 3] = 1.0 - phiA[t];
      ps[2, i, t, 1] = phiB[t] * psiBA[t];
      ps[2, i, t, 2] = phiB[t] * (1 - psiBA[t]);
      ps[2, i, t, 3] = 1.0 - phiB[t];
      ps[3, i, t, 1] = 0.0;
      ps[3, i, t, 2] = 0.0;
      ps[3, i, t, 3] = 1.0;

      // Define probabilities of O(t) given S(t)
      po[1, i, t, 1] = pA[t];
      po[1, i, t, 2] = 0.0;
      po[1, i, t, 3] = 1.0 - pA[t];
      po[2, i, t, 1] = 0.0;
      po[2, i, t, 2] = pB[t];
      po[2, i, t, 3] = 1.0 - pB[t];
      po[3, i, t, 1] = 0.0;
      po[3, i, t, 2] = 0.0;
      po[3, i, t, 3] = 1.0;
    }
  }
}

model {
  real acc[3];
  vector[3] gamma[n_occasions];

  // Priors
  // Uniform priors are implicitly defined.
  //  mean_phi ~ uniform(0, 1);
  //  mean_psi ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);

  // Likelihood
  // Forward algorithm derived from Stan Modeling Language
  // User's Guide and Reference Manual
  for (i in 1:nind) {
    if (first[i] > 0) {
      for (k in 1:3)
        gamma[first[i], k] = (k == y[i, first[i]]);

      for (t in (first[i] + 1):n_occasions) {
        for (k in 1:3) {
          for (j in 1:3)
            acc[j] = gamma[t - 1, j] * ps[j, i, t - 1, k]
                    * po[k, i, t - 1, y[i, t]];
          gamma[t, k] = sum(acc);
        }
      }
      target += log(sum(gamma[n_occasions]));
    }
  }
}
