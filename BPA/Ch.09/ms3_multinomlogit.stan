// -------------------------------------------------
// States (S):
// 1 alive at A
// 2 alive at B
// 3 alive at C
// 4 dead
// Observations (O):
// 1 seen at A
// 2 seen at B
// 3 seen at C
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

  /**
   * Return a simplex such as follows (thanks to Bob Carpenter):
   * p[1] <- exp(lp[1]) / (1.0 + exp(lp[1]) + exp(lp[2]));
   * p[2] <- exp(lp[2]) / (1.0 + exp(lp[1]) + exp(lp[2]));
   * p[3] <- 1.0 - p[1] - p[2];
   *
   * @param lp   N-dimension vector
   * @return (N+1)-simplex of given vector and 0
   */
  vector softmax_0(vector lp) {
    vector[num_elements(lp) + 1] lp_temp;

    lp_temp[1:num_elements(lp)] = lp;
    lp_temp[num_elements(lp) + 1] = 0;
    return softmax(lp_temp);
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
  real<lower=0,upper=1> phiA; // Survival probability at site A
  real<lower=0,upper=1> phiB; // Survival probability at site B
  real<lower=0,upper=1> phiC; // Survival probability at site C
  real<lower=0,upper=1> pA;   // Recapture probability at site A
  real<lower=0,upper=1> pB;   // Recapture probability at site B
  real<lower=0,upper=1> pC;   // Recapture probability at site C
  vector[2] lpsiA;            // Logit of movement probability from site A
  vector[2] lpsiB;            // Logit of movement probability from site B
  vector[2] lpsiC;            // Logit of movement probability from site C
}

transformed parameters {
  simplex[3] psiA; // Movement probability from site A
  simplex[3] psiB; // Movement probability from site B
  simplex[3] psiC; // Movement probability from site C
  simplex[4] ps[4, nind, n_occ_minus_1];
  simplex[4] po[4, nind, n_occ_minus_1];

  // Constrain the transitions such that their sum is < 1
  psiA = softmax_0(lpsiA);
  psiB = softmax_0(lpsiB);
  psiC = softmax_0(lpsiC);

  // Define state-transition and observation matrices
  for (i in 1:nind) {
    // Define probabilities of state S(t+1) given S(t)
    for (t in 1:n_occ_minus_1) {
      ps[1, i, t, 1] = phiA * psiA[1];
      ps[1, i, t, 2] = phiA * psiA[2];
      ps[1, i, t, 3] = phiA * psiA[3];
      ps[1, i, t, 4] = 1.0 - phiA;
      ps[2, i, t, 1] = phiB * psiB[1];
      ps[2, i, t, 2] = phiB * psiB[2];
      ps[2, i, t, 3] = phiB * psiB[3];
      ps[2, i, t, 4] = 1.0 - phiB;
      ps[3, i, t, 1] = phiC * psiC[1];
      ps[3, i, t, 2] = phiC * psiC[2];
      ps[3, i, t, 3] = phiC * psiC[3];
      ps[3, i, t, 4] = 1.0 - phiC;
      ps[4, i, t, 1] = 0.0;
      ps[4, i, t, 2] = 0.0;
      ps[4, i, t, 3] = 0.0;
      ps[4, i, t, 4] = 1.0;

      // Define probabilities of O(t) given S(t)
      po[1, i, t, 1] = pA;
      po[1, i, t, 2] = 0.0;
      po[1, i, t, 3] = 0.0;
      po[1, i, t, 4] = 1.0 - pA;
      po[2, i, t, 1] = 0.0;
      po[2, i, t, 2] = pB;
      po[2, i, t, 3] = 0.0;
      po[2, i, t, 4] = 1.0 - pB;
      po[3, i, t, 1] = 0.0;
      po[3, i, t, 2] = 0.0;
      po[3, i, t, 3] = pC;
      po[3, i, t, 4] = 1.0 - pC;
      po[4, i, t, 1] = 0.0;
      po[4, i, t, 2] = 0.0;
      po[4, i, t, 3] = 0.0;
      po[4, i, t, 4] = 1.0;
      }
   }
}

model {
  real acc[4];
  vector[4] gamma[n_occasions];

  // Priors
  // Survival and recapture: uniform
  // Uniform priors are implicitly defined.
  //  phiA ~ uniform(0, 1);
  //  phiB ~ uniform(0, 1);
  //  phiC ~ uniform(0, 1);
  //  pA ~ uniform(0, 1);
  //  pB ~ uniform(0, 1);
  //  pC ~ uniform(0, 1);

  // Normal priors on logit of all but one transition probs
  lpsiA ~ normal(0, sqrt(1000));
  lpsiB ~ normal(0, sqrt(1000));
  lpsiC ~ normal(0, sqrt(1000));

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
