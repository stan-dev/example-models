data {
  int<lower=0> n_occasions;
  int<lower=0> marr[n_occasions-1,n_occasions];
}

transformed data {
  int r[n_occasions-1];

  // Calculate the number of birds released each year
  for (t in 1:(n_occasions - 1))
     r[t] <- sum(marr[t]);
}

parameters {
  vector<lower=0,upper=1>[n_occasions - 1] phi; // Survival
  vector<lower=0,upper=1>[n_occasions - 1] p;   // Recapture
}

transformed parameters {
  vector<lower=0,upper=1>[n_occasions - 1] q;
  simplex[n_occasions] pr[n_occasions - 1];

  q <- 1.0 - p;             // Probability of non-recapture

  // Define the cell probabilities of the m-array
  // Main diagonal
  for (t in 1:(n_occasions - 1)) {
    pr[t][t] <- phi[t] * p[t];
    // Above main diagonal
      for (j in (t + 1):(n_occasions - 1)) {
        pr[t][j] <- prod(segment(phi, t, j - t + 1)) *
                    prod(segment(q, t, j - t)) *
                    p[j];
      } //j
    // Below main diagonal
    for (j in 1:(t - 1))
      pr[t][j] <- 0;
  } //t

  // Last column: probability of non-recapture
  for (t in 1:(n_occasions - 1))
    pr[t][n_occasions] <- 1 - sum(head(pr[t], n_occasions - 1));
}

model {
  // Priors
  phi ~ uniform(0, 1);
  p ~ uniform(0, 1);

  // Define the multinomial likelihood
  for (t in 1:(n_occasions - 1)) {
    increment_log_prob(multinomial_log(marr[t], pr[t]));
  }
}

generated quantities {
  real fit;
  real fit_new;
  matrix[n_occasions - 1, n_occasions] E_org;
  matrix[n_occasions - 1, n_occasions] E_new;
  vector[n_occasions] expmarr[n_occasions - 1];
  int<lower=0> marr_new[n_occasions - 1, n_occasions];

  // Assess model fit using Freeman-Tukey statistic
  // Compute fit statistics for observed data
  for (t in 1:(n_occasions - 1)){
    expmarr[t] <- r[t] * pr[t];
    for (j in 1:n_occasions){
      E_org[t, j] <- square((sqrt(marr[t][j]) - sqrt(expmarr[t][j])));
    } //j
  } //t

  // Generate replicate data and compute fit stats from them
  for (t in 1:(n_occasions - 1)) {
    marr_new[t] <- multinomial_rng(pr[t], r[t]);
    for (j in 1:n_occasions) {
      E_new[t, j] <- square((sqrt(marr_new[t][j]) - sqrt(expmarr[t][j])));
    }
  }
  fit <- sum(E_org);
  fit_new <- sum(E_new);
}
