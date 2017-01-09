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
  int<lower=0,upper=n_occasions> first[nind];

  for (i in 1:nind)
    first[i] = first_capture(y[i]);
}

parameters {
  real<lower=0,upper=1> phiA;   // Mean survival in A
  real<lower=0,upper=1> phiB;   // Mean survival in B
  real<lower=0,upper=1> psiAB;  // Mean movement from A to B
  real<lower=0,upper=1> psiBA;  // Mean movement from B to A
  real<lower=0,upper=1> pA;     // Mean recapture in A
  real<lower=0,upper=1> pB;     // Mean recapture in B
}

transformed parameters {
  simplex[3] ps[3];
  simplex[3] po[3];

  // Define state-transition and observation matrices
  // Define probabilities of state S(t+1) given S(t)
  ps[1, 1] = phiA * (1.0 - psiAB);
  ps[1, 2] = phiA * psiAB;
  ps[1, 3] = 1.0 - phiA;
  ps[2, 1] = phiB * psiBA;
  ps[2, 2] = phiB * (1 - psiBA);
  ps[2, 3] = 1.0 - phiB;
  ps[3, 1] = 0.0;
  ps[3, 2] = 0.0;
  ps[3, 3] = 1.0;

  // Define probabilities of O(t) given S(t)
  po[1, 1] = pA;
  po[1, 2] = 0.0;
  po[1, 3] = 1.0 - pA;
  po[2, 1] = 0.0;
  po[2, 2] = pB;
  po[2, 3] = 1.0 - pB;
  po[3, 1] = 0.0;
  po[3, 2] = 0.0;
  po[3, 3] = 1.0;
}

model {
  real acc[3];
  vector[3] gamma[n_occasions];

  // Priors
  // Uniform priors are implicitly defined.
  //  phiA ~ uniform(0, 1);
  //  phiB ~ uniform(0, 1);
  //  psiAB ~ uniform(0, 1);
  //  psiBA ~ uniform(0, 1);
  //  pA ~ uniform(0, 1);
  //  pB ~ uniform(0, 1);

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
            acc[j] = gamma[t - 1, j] * ps[j, k]
                    * po[k, y[i, t]];
          gamma[t, k] = sum(acc);
        }
      }
      target += log(sum(gamma[n_occasions]));
    }
  }
}
