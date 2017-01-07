data {
  int n_occasions;
  int marr[n_occasions, n_occasions + 1];
}

parameters {
  real<lower=0,upper=1> mean_s;           // Mean survival
  real<lower=0,upper=1> mean_r;           // Mean recovery
}

transformed parameters {
  vector[n_occasions] s;
  vector[n_occasions] r;
  simplex[n_occasions + 1] pr[n_occasions];

  // Constraints
  s = rep_vector(mean_s, n_occasions);
  r = rep_vector(mean_r, n_occasions);

  // Define the cell probabilities of the m-array
  // Main diagonal
  for (t in 1:n_occasions) {
    pr[t, t] = (1 - s[t]) * r[t];

    // Above main diagonal
    for (j in (t + 1):n_occasions)
      pr[t, j] = prod(s[t:(j - 1)]) * (1 - s[j]) * r[j];

    // Below main diagonal
    for (j in 1:(t - 1))
      pr[t, j] = 0;
  }

  // Last column: probability of non-recovery
  for (t in 1:n_occasions)
    pr[t, n_occasions + 1] = 1 - sum(pr[t, :n_occasions]);
}

model {
  // Priors
  // Uniform priors are implicitly defined.
  //  mean_s ~ uniform(0, 1);
  //  mean_r ~ uniform(0, 1);

  // Define the multinomial likelihood
  for (t in 1:n_occasions)
    marr[t] ~ multinomial(pr[t]);
}
