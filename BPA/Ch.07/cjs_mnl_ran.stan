data {
  int<lower=0> n_occasions;
  int<lower=0> marr[n_occasions - 1, n_occasions];
}

transformed data {
  int r[n_occasions];

  for (t in 1:n_occasions - 1)
    r[t] <- sum(marr[t]);
}

parameters {
  real<lower=0,upper=1> mean_phi;       // Mean survival
  real<lower=0,upper=1> mean_p;         // Mean recapture
  real<lower=0> sigma;
  real epsilon[n_occasions - 1];
}

transformed parameters {
  real mu;
  vector<lower=0,upper=1>[n_occasions - 1] phi;
  vector<lower=0,upper=1>[n_occasions - 1] p;
  vector<lower=0,upper=1>[n_occasions - 1] q;
  simplex[n_occasions] pr[n_occasions - 1];

  mu <- logit(mean_phi);

  // Constraints
  for (t in 1:(n_occasions - 1)) {
    phi[t] <- inv_logit(mu + epsilon[t]);
    p[t] <- mean_p;
  }

  q <- 1.0 - p;             // Probability of non-recapture

  // Define the cell probabilities of the m-arrays
  // Main diagonal
  for (t in 1:(n_occasions - 1)) {
    pr[t][t] <- phi[t] * p[t];

    // Above main diagonal
    for (j in (t + 1):(n_occasions - 1))
       pr[t][j] <- prod(segment(phi, t, j - t + 1)) *
                   prod(segment(q, t, j - t)) *
                   p[j];

    // Below main diagonal
    for (j in 1:(t - 1))
      pr[t][j] <- 0;
  } // t

  // Last column: probability of non-recapture
  for (t in 1:(n_occasions - 1))
    pr[t][n_occasions] <- 1 - sum(head(pr[t], n_occasions - 1));
}

model {
  // Priors
  epsilon ~ normal(0, sigma);
  mean_phi ~ uniform(0, 1);
  sigma ~ uniform(0, 5);
  mean_p ~ uniform(0, 1);

  // Define the multinomial likelihood
  for (t in 1:(n_occasions - 1))
    increment_log_prob(multinomial_log(marr[t], pr[t]));
}

generated quantities {
  real<lower=0> sigma2;
  real<lower=0> sigma2_real;
  vector[n_occasions] expmarr[n_occasions - 1];
  int marr_new[n_occasions - 1, n_occasions];
  matrix[n_occasions - 1, n_occasions] E_org;
  matrix[n_occasions - 1, n_occasions] E_new;
  real fit;
  real fit_new;

  sigma2 <- square(sigma);

  // Temporal variance on real scale
  sigma2_real <- sigma2 * square(mean_phi) * square(1 - mean_phi);

  // Assess model fit using Freeman-Tukey statistic
  // Compute fit statistics for observed data
  for (t in 1:(n_occasions - 1)) {
    expmarr[t] <- r[t] * pr[t];
    for (j in 1:n_occasions)
      E_org[t, j] <- square((sqrt(marr[t, j]) - sqrt(expmarr[t][j])));
   }

  // Generate replicate data and compute fit stats from them
  for (t in 1:(n_occasions - 1)) {
    marr_new[t] <- multinomial_rng(pr[t], r[t]);
    for (j in 1:n_occasions)
      E_new[t, j] <- square((sqrt(marr_new[t, j]) - sqrt(expmarr[t][j])));
  }
  fit <- sum(E_org);
  fit_new <- sum(E_new);
}
