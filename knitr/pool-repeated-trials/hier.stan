data {
  int<lower=0> N;              // items
  int<lower=0> K[N];           // initial trials
  int<lower=0> y[N];           // initial successes

  int<lower=0> K_new[N];       // new trials
  int<lower=0> y_new[N];       // new successes
}
transformed data {
  real min_y;   // minimum successes
  real max_y;   // maximum successes
  real mean_y;  // sample mean successes
  real sd_y;    // sample std dev successes

  min_y <- min(y);
  max_y <- max(y);
  mean_y <- mean(to_vector(y));
  sd_y <- sd(to_vector(y));
}
parameters {
  real<lower=0, upper=1> phi;         // population chance of success
  real<lower=1> kappa;                // population concentration
  vector<lower=0, upper=1>[N] theta;  // chance of success 
}
model {
  kappa ~ pareto(1, 1.5);                        // hyperprior
  theta ~ beta(phi * kappa, (1 - phi) * kappa);  // prior
  y ~ binomial(K, theta);                        // likelihood
}
generated quantities {
  real log_p_new;     // posterior predictive log density remaining trials

  int<lower=0> z[N];  // posterior prediction remaining trials

  int<lower=0, upper=1> some_ability_gt_350;  // Pr[some theta > 0.35]
  int<lower=0, upper=1> avg_gt_400[N];        // Pr[season avg of n] >= 0.400
  int<lower=0, upper=1> ability_gt_400[N];    // Pr[chance-of-success of n] >= 0.400

  int<lower=1, upper=N> rnk[N];      // rank of player n
  int<lower=0, upper=1> is_best[N];  // Pr[player n highest chance of success]

  int<lower=0> y_rep[N];      // replications for existing items
  int<lower=0> y_pop_rep[N];  // replications for simulated items

  real min_y_rep;   // posterior predictive min replicated successes
  real max_y_rep;   // posterior predictive max replicated successes
  real mean_y_rep;  // posterior predictive sample mean replicated successes
  real sd_y_rep;    // posterior predictive sample std dev replicated successes

  int p_min;                                     // posterior predictive p-values
  int p_max;
  int p_mean;
  int p_sd;

  log_p_new <- 0;
  for (n in 1:N)
    log_p_new <- log_p_new + binomial_log(y_new[n], K_new[n], theta[n]);
  
  for (n in 1:N)
    z[n] <- binomial_rng(K_new[n], theta[n]);

  some_ability_gt_350 <- (max(theta) > 0.35);
  for (n in 1:N)
    avg_gt_400[n] <- (((y[n] + z[n]) / (0.0 + K [n] + K_new[n])) > 0.400);
  for (n in 1:N)
    ability_gt_400[n] <- (theta[n] > 0.400);

  {
    int dsc[N];
    dsc <- sort_indices_desc(theta);
    for (n in 1:N)
      rnk[dsc[n]] <- n;
  }
  for (n in 1:N)
    is_best[n] <- (rnk[n] == 1);

  for (n in 1:N)
    y_rep[n] <- binomial_rng(K[n], theta[n]);
  for (n in 1:N)
    y_pop_rep[n] <- binomial_rng(K[n], 
                                 beta_rng(phi * kappa,
                                          (1 - phi) * kappa));

  min_y_rep <- min(y_rep);
  max_y_rep <- max(y_rep);
  mean_y_rep <- mean(to_vector(y_rep));
  sd_y_rep <- sd(to_vector(y_rep));

  p_min <- (min_y_rep >= min_y);
  p_max <- (max_y_rep >= max_y);
  p_mean <- (mean_y_rep >= mean_y);
  p_sd <- (sd_y_rep >= sd_y);
}
