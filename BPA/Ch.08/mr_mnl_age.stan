data {
  int n_occasions;
  int marr_j[n_occasions, n_occasions + 1];
  int marr_a[n_occasions, n_occasions + 1];
}

parameters {
  real<lower=0,upper=1> mean_sj;         // Mean juv. survival
  real<lower=0,upper=1> mean_sa;         // Mean ad. survival
  real<lower=0,upper=1> mean_rj;         // Mean juv. recovery
  real<lower=0,upper=1> mean_ra;         // Mean ad. recovery
}

transformed parameters {
  vector[n_occasions] sj;
  vector[n_occasions] sa;
  vector[n_occasions] rj;
  vector[n_occasions] ra;
  simplex[n_occasions + 1] pr_a[n_occasions];
  simplex[n_occasions + 1] pr_j[n_occasions];

  // Constraints
  for (t in 1:n_occasions) {
    sj[t] <- mean_sj;
    sa[t] <- mean_sa;
    rj[t] <- mean_rj;
    ra[t] <- mean_ra;
  }

  // Define the cell probabilities of the juvenile m-array
  // Main diagonal
  for (t in 1:n_occasions) {
    pr_j[t, t] <- (1 - sj[t]) * rj[t];

    // Further above main diagonal
    for (j in (t + 2):n_occasions)
      pr_j[t, j] <- sj[t] * prod(segment(sa, t + 1, j - t - 1)) *
                    (1 - sa[j]) * ra[j];

    // Below main diagonal
    for (j in 1:(t - 1))
      pr_j[t, j] <- 0;
  } //t

  for (t in 1:(n_occasions - 1))
    // One above main diagonal
    pr_j[t, t + 1] <- sj[t] * (1 - sa[t + 1]) * ra[t + 1];

  // Last column: probability of non-recovery
  for (t in 1:n_occasions)
    pr_j[t, n_occasions + 1] <- 1 - sum(head(pr_j[t], n_occasions));

  // Define the cell probabilities of the adult m-array
  // Main diagonal
  for (t in 1:n_occasions) {
    pr_a[t,t] <- (1 - sa[t]) * ra[t];
    // Above main diagonal
    for (j in (t + 1):n_occasions)
      pr_a[t, j] <- prod(segment(sa, t, j - t)) * (1 - sa[j]) * ra[j];

    // Below main diagonal
    for (j in 1:(t - 1))
      pr_a[t, j] <- 0;
  } //t

  // Last column: probability of non-recovery
  for (t in 1:n_occasions)
   pr_a[t, n_occasions + 1] <- 1 - sum(head(pr_a[t], n_occasions));
}

model {
  // Priors
  mean_sj ~ uniform(0, 1);
  mean_sa ~ uniform(0, 1);
  mean_rj ~ uniform(0, 1);
  mean_ra ~ uniform(0, 1);

  // Define the multinomial likelihood
  for (t in 1:n_occasions) {
    increment_log_prob(multinomial_log(marr_j[t], pr_j[t]));
    increment_log_prob(multinomial_log(marr_a[t], pr_a[t]));
  }
}
