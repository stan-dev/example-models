// ---------------------------------
// States (S):
// 1 alive and present
// 2 alive and absent
// 3 dead
// Observations (O):
// 1 seen
// 2 not seen
// ---------------------------------

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
      if (y_i[k] == 1)
        return k;
    return 0;
  }
}

data {
  int<lower=0> nind;
  int<lower=0> n_occasions;
  int<lower=1,upper=2> y[nind, n_occasions];
}

transformed data {
  int n_occ_minus_1 = n_occasions - 1;
  int<lower=0,upper=n_occasions> first[nind];

  for (i in 1:nind)
    first[i] = first_capture(y[i]);
}

parameters {
  real<lower=0,upper=1> mean_phi;   // Mean state-spec. survival
  real<lower=0,upper=1> mean_psiIO; // Mean temp. emigration
  real<lower=0,upper=1> mean_psiOI; // Mean temp. immigration
  real<lower=0,upper=1> mean_p;     // Mean state-spec. recapture
}

transformed parameters {
  vector<lower=0,upper=1>[n_occ_minus_1] phi;   // Survival probability
  vector<lower=0,upper=1>[n_occ_minus_1] psiIO; // Probability to emigrate
  vector<lower=0,upper=1>[n_occ_minus_1] psiOI; // Probability to immigrate
  vector<lower=0,upper=1>[n_occ_minus_1] p;     // Recapture probability
  simplex[3] ps[3, nind, n_occ_minus_1];
  simplex[2] po[3, nind, n_occ_minus_1];

  // Constraints
  for (t in 1:n_occ_minus_1) {
    phi[t] = mean_phi;
    psiIO[t] = mean_psiIO;
    psiOI[t] = mean_psiOI;
    p[t] = mean_p;
  }

  // Define state-transition and observation matrices
  for (i in 1:nind) {
    // Define probabilities of state S(t+1) given S(t)
    for (t in 1:n_occ_minus_1) {
      ps[1, i, t, 1] = phi[t] * (1.0 - psiIO[t]);
      ps[1, i, t, 2] = phi[t] * psiIO[t];
      ps[1, i, t, 3] = 1.0 - phi[t];
      ps[2, i, t, 1] = phi[t] * psiOI[t];
      ps[2, i, t, 2] = phi[t] * (1.0 - psiOI[t]);
      ps[2, i, t, 3] = 1.0 - phi[t];
      ps[3, i, t, 1] = 0.0;
      ps[3, i, t, 2] = 0.0;
      ps[3, i, t, 3] = 1.0;

      // Define probabilities of O(t) given S(t)
      po[1, i, t, 1] = p[t];
      po[1, i, t, 2] = 1.0 - p[t];
      po[2, i, t, 1] = 0.0;
      po[2, i, t, 2] = 1.0;
      po[3, i, t, 1] = 0.0;
      po[3, i, t, 2] = 1.0;
      }
   }
}

model {
  real acc[3];
  vector[3] gamma[n_occasions];

  // Priors
  // Uniform priors are implicitly defined.
  //  mean_phi ~ uniform(0, 1);
  //  mean_psiIO ~ uniform(0, 1);
  //  mean_psiOI ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);

  // Likelihood
  // Forward algorithm derived from Stan Modeling Language
  // User's Guide and Reference Manual.
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
