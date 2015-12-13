data {
  int<lower=0> n_occasions;
  int<lower=0> marr_j[n_occasions - 1, n_occasions];
  int<lower=0> marr_a[n_occasions - 1, n_occasions];
}

parameters {
  real<lower=0,upper=1> mean_phijuv;     // Mean juv. survival
  real<lower=0,upper=1> mean_phiad;      // Mean ad. survival
  real<lower=0,upper=1> mean_p;          // Mean recapture
}

transformed parameters {
  vector<lower=0,upper=1>[n_occasions - 1] phi_juv;
  vector<lower=0,upper=1>[n_occasions - 1] phi_ad;
  vector<lower=0,upper=1>[n_occasions - 1] p;
  vector<lower=0,upper=1>[n_occasions - 1] q;
  simplex[n_occasions] pr_j[n_occasions - 1];
  simplex[n_occasions] pr_a[n_occasions - 1];

  // Constraints
  for (t in 1:(n_occasions - 1)) {
    phi_juv[t] <- mean_phijuv;
    phi_ad[t] <- mean_phiad;
    p[t] <- mean_p;
  }

  q <- 1.0 - p;             // Probability of non-recapture

  // Define the cell probabilities of the m-arrays
  // Main diagonal
  for (t in 1:(n_occasions - 1)) {
    pr_j[t][t] <- phi_juv[t] * p[t];
    pr_a[t][t] <- phi_ad[t] * p[t];

    // Above main diagonal
    for (j in (t + 1):(n_occasions - 1)) {
       pr_j[t][j] <- phi_juv[t] *
                     prod(segment(phi_ad, t + 1, j - t)) *
                     prod(segment(q, t, j - t)) *
                     p[j];
       pr_a[t][j] <- prod(segment(phi_ad, t, j - t + 1)) *
                     prod(segment(q, t, j - t)) *
                     p[j];
    } // j
    // Below main diagonal
    for (j in 1:(t - 1)) {
      pr_j[t][j] <- 0;
      pr_a[t][j] <- 0;
    } //j
  } // t

  // Last column: probability of non-recapture
  for (t in 1:(n_occasions - 1)) {
    pr_j[t][n_occasions] <- 1 - sum(head(pr_j[t], n_occasions - 1));
    pr_a[t][n_occasions] <- 1 - sum(head(pr_a[t], n_occasions - 1));
  } // t
}

model {
  // Priors
  mean_phijuv ~ uniform(0, 1);
  mean_phiad ~ uniform(0, 1);
  mean_p ~ uniform(0, 1);

  // Define the multinomial likelihood
  for (t in 1:(n_occasions - 1)) {
    increment_log_prob(multinomial_log(marr_j[t], pr_j[t]));
    increment_log_prob(multinomial_log(marr_a[t], pr_a[t]));
   }
}
