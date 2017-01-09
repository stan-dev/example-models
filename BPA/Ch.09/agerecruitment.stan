// -------------------------------------------------
// States (S):
// 1 juvenile
// 2 not yet breeding at age 1 year
// 3 not yet breeding at age 2 years
// 4 breeder
// 5 dead
// Observations (O):
// 1 seen as juvenile
// 2 seen as not yet breeding
// 3 seen breeding
// 4 not seen
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
      if (y_i[k] != 4)
        return k;
    return 0;
  }
}

data {
  int<lower=0> nind;
  int<lower=0> n_occasions;
  int<lower=1,upper=4> y[nind, n_occasions];
}

transformed data {
  int n_occ_minus_1 = n_occasions - 1;
  int<lower=0,upper=n_occasions> first[nind];

  for (i in 1:nind)
    first[i] = first_capture(y[i]);
}

parameters {
  real<lower=0,upper=1> mean_phi1;    // Mean 1y survival
  real<lower=0,upper=1> mean_phi2;    // Mean 2y survival
  real<lower=0,upper=1> mean_phiad;   // Mean ad survival
  real<lower=0,upper=1> mean_alpha1;  // Mean 1y breeding prob.
  real<lower=0,upper=1> mean_alpha2;  // Mean 2y breeding prob.
  real<lower=0,upper=1> mean_pNB;     // Mean recapture non-breeders
  real<lower=0,upper=1> mean_pB;      // Mean recapture breeders
}

transformed parameters {
  vector<lower=0,upper=1>[n_occ_minus_1] phi_1;   // First year survival prob.
  vector<lower=0,upper=1>[n_occ_minus_1] phi_2;   // Second year survival prob.
  vector<lower=0,upper=1>[n_occ_minus_1] phi_ad;  // Adult survival prob.
  vector<lower=0,upper=1>[n_occ_minus_1] alpha_1; // Prob. to start breeding when 1 yr old
  vector<lower=0,upper=1>[n_occ_minus_1] alpha_2; // Prob. to start breeding when 2 yr old
  vector<lower=0,upper=1>[n_occ_minus_1] p_NB;    // Recapture prob. of non-breeders
  vector<lower=0,upper=1>[n_occ_minus_1] p_B;     // Recapture prob. of breeders
  simplex[5] ps[5, nind, n_occ_minus_1];
  simplex[4] po[5, nind, n_occ_minus_1];

  // Constraints
  for (t in 1:n_occ_minus_1) {
    phi_1[t] = mean_phi1;
    phi_2[t] = mean_phi2;
    phi_ad[t] = mean_phiad;
    alpha_1[t] = mean_alpha1;
    alpha_2[t] = mean_alpha2;
    p_NB[t] = mean_pNB;
    p_B[t] = mean_pB;
  }

  // Define state-transition and observation matrices
  for (i in 1:nind) {
    // Define probabilities of state S(t+1) given S(t)
    for (t in 1:n_occ_minus_1) {
      ps[1, i, t, 1] = 0.0;
      ps[1, i, t, 2] = phi_1[t] * (1 - alpha_1[t]);
      ps[1, i, t, 3] = 0.0;
      ps[1, i, t, 4] = phi_1[t] * alpha_1[t];
      ps[1, i, t, 5] = 1.0 - phi_1[t];
      ps[2, i, t, 1] = 0.0;
      ps[2, i, t, 2] = 0.0;
      ps[2, i, t, 3] = phi_2[t] * (1.0 - alpha_2[t]);
      ps[2, i, t, 4] = phi_2[t] * alpha_2[t];
      ps[2, i, t, 5] = 1.0 - phi_2[t];
      ps[3, i, t, 1] = 0.0;
      ps[3, i, t, 2] = 0.0;
      ps[3, i, t, 3] = 0.0;
      ps[3, i, t, 4] = phi_ad[t];
      ps[3, i, t, 5] = 1 - phi_ad[t];
      ps[4, i, t, 1] = 0.0;
      ps[4, i, t, 2] = 0.0;
      ps[4, i, t, 3] = 0.0;
      ps[4, i, t, 4] = phi_ad[t];
      ps[4, i, t, 5] = 1.0 - phi_ad[t];
      ps[5, i, t, 1] = 0.0;
      ps[5, i, t, 2] = 0.0;
      ps[5, i, t, 3] = 0.0;
      ps[5, i, t, 4] = 0.0;
      ps[5, i, t, 5] = 1.0;

      // Define probabilities of O(t) given S(t)
      po[1, i, t, 1] = 0.0;
      po[1, i, t, 2] = 0.0;
      po[1, i, t, 3] = 0.0;
      po[1, i, t, 4] = 1.0;
      po[2, i, t, 1] = 0.0;
      po[2, i, t, 2] = p_NB[t];
      po[2, i, t, 3] = 0.0;
      po[2, i, t, 4] = 1.0 - p_NB[t];
      po[3, i, t, 1] = 0.0;
      po[3, i, t, 2] = p_NB[t];
      po[3, i, t, 3] = 0.0;
      po[3, i, t, 4] = 1.0 - p_NB[t];
      po[4, i, t, 1] = 0.0;
      po[4, i, t, 2] = 0.0;
      po[4, i, t, 3] = p_B[t];
      po[4, i, t, 4] = 1.0 - p_B[t];
      po[5, i, t, 1] = 0.0;
      po[5, i, t, 2] = 0.0;
      po[5, i, t, 3] = 0.0;
      po[5, i, t, 4] = 1.0;
    }
  }
}

model {
  real acc[5];
  vector[5] gamma[n_occasions];

  // Priors
  // Uniform priors are implicitly defined.
  //  mean_phi1 ~ uniform(0, 1);
  //  mean_phi2 ~ uniform(0, 1);
  //  mean_phiad ~ uniform(0, 1);
  //  mean_alpha1 ~ uniform(0, 1);
  //  mean_alpha2 ~ uniform(0, 1);
  //  mean_pNB ~ uniform(0, 1);
  //  mean_pB ~ uniform(0, 1);

  // Likelihood
  // Forward algorithm derived from Stan Modeling Language
  // User's Guide and Reference Manual
  for (i in 1:nind) {
    if (first[i] > 0) {
      for (k in 1:5)
        gamma[first[i], k] = (k == y[i, first[i]]);

      for (t in (first[i] + 1):n_occasions) {
        for (k in 1:5) {
          for (j in 1:5)
            acc[j] = gamma[t - 1, j] * ps[j, i, t - 1, k]
                    * po[k, i, t - 1, y[i, t]];
          gamma[t, k] = sum(acc);
        }
      }
      target += log(sum(gamma[n_occasions]));
    }
  }
}
