data {
  int<lower=0> n_occasions;     // Number of capture occasions
  int<lower=0> marr_j[n_occasions - 1, n_occasions];    // Juv. m-array
  int<lower=0> marr_a[n_occasions - 1, n_occasions];    // Ad. m-array
}

transformed data {
  // Compoud declaration is enabled in Stan 2.13
  int n_occ_minus_1 = n_occasions - 1;
  //  int n_occ_minus_1;

  //  n_occ_minus_1 = n_occasions - 1;
}

parameters {
  real<lower=0,upper=1> mean_phijuv;     // Mean juv. survival
  real<lower=0,upper=1> mean_phiad;      // Mean ad. survival
  real<lower=0,upper=1> mean_p;          // Mean recapture
}

transformed parameters {
  vector<lower=0,upper=1>[n_occ_minus_1] phi_juv;
  vector<lower=0,upper=1>[n_occ_minus_1] phi_ad;
  vector<lower=0,upper=1>[n_occ_minus_1] p;
  vector<lower=0,upper=1>[n_occ_minus_1] q;
  simplex[n_occasions] pr_j[n_occ_minus_1];
  simplex[n_occasions] pr_a[n_occ_minus_1];

  // Constraints
  phi_juv = rep_vector(mean_phijuv, n_occ_minus_1);
  phi_ad = rep_vector(mean_phiad, n_occ_minus_1);
  p = rep_vector(mean_p, n_occ_minus_1);

  q = 1.0 - p;             // Probability of non-recapture

  // Define the cell probabilities of the m-arrays
  // Main diagonal
  for (t in 1:n_occ_minus_1) {
    pr_j[t, t] = phi_juv[t] * p[t];
    pr_a[t, t] = phi_ad[t] * p[t];

    // Above main diagonal
    for (j in (t + 1):n_occ_minus_1) {
      pr_j[t, j] = phi_juv[t] * prod(phi_ad[(t + 1):j])
                  * prod(q[t:(j - 1)]) * p[j];
      pr_a[t, j] = prod(phi_ad[t:j]) * prod(q[t:(j - 1)]) * p[j];
    }

    // Below main diagonal
    pr_j[t, :(t - 1)] = rep_vector(0, t - 1);
    pr_a[t, :(t - 1)] = rep_vector(0, t - 1);
  }

  // Last column: probability of non-recapture
  for (t in 1:n_occ_minus_1) {
    pr_j[t, n_occasions] = 1 - sum(pr_j[t, :n_occ_minus_1]);
    pr_a[t, n_occasions] = 1 - sum(pr_a[t, :n_occ_minus_1]);
  }
}

model {
  // Priors
  // Uniform priors are implicitly defined.
  //  mean_phijuv ~ uniform(0, 1);
  //  mean_phiad ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);

  // Define the multinomial likelihood
  for (t in 1:n_occ_minus_1) {
    marr_j[t] ~ multinomial(pr_j[t]);
    marr_a[t] ~ multinomial(pr_a[t]);
   }
}
