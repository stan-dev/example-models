data {
  int<lower=0> N; // items
  array[N] int<lower=0> K; // initial trials
  array[N] int<lower=0> y; // initial successes
  
  array[N] int<lower=0> K_new; // new trials
  array[N] int<lower=0> y_new; // new successes
}
transformed data {
  real min_y; // minimum successes
  real max_y; // maximum successes
  real mean_y; // sample mean successes
  real sd_y; // sample std dev successes
  
  min_y = min(y);
  max_y = max(y);
  mean_y = mean(to_vector(y));
  sd_y = sd(to_vector(y));
}
parameters {
  real<lower=0, upper=1> phi; // chance of success (pooled)
}
model {
  y ~ binomial(K, phi); // likelihood
}
generated quantities {
  // per-player chance-of-success
  vector<lower=0, upper=1>[N] theta = rep_vector(phi, N);

  // posterior predictive log density remaining trials (pooled)
  real log_p_new = 0; 
  for (n in 1 : N) {
    log_p_new = log_p_new + binomial_lpmf(y_new[n] | K_new[n], theta[n]);
  }
  
  // posterior predictive remaining trials
  array[N] int<lower=0> z = binomial_rng(K_new, theta); 

  // Pr[theta > 0.35] (pooled)
  int<lower=0, upper=1> some_ability_gt_350 = phi > 0.35; 

  // Pr[some player avgerage/ability > 400]
  array[N] int<lower=0, upper=1> avg_gt_400;
  for (n in 1 : N) {
    avg_gt_400[n] = ((y[n] + z[n]) / (0.0 + K[n] + K_new[n])) > 0.400;
  }
  array[N] int<lower=0, upper=1> ability_gt_400;
  for (n in 1 : N) {
    ability_gt_400[n] = theta[n] > 0.400;
  }

  // replications for existing items  
  array[N] int<lower=0> y_rep = binomial_rng(K, theta);
  
  // posterior predictive min replicated successes
  real<lower=0> min_y_rep = min(y_rep);
  int<lower=0, upper=1> p_min = min_y_rep >= min_y;

  // posterior predictive max replicated successes
  real<lower=0> max_y_rep = max(y_rep); 
  int<lower=0, upper=1> p_max = max_y_rep >= max_y;

  // posterior predictive sample mean replicated successes
  real<lower=0> mean_y_rep = mean(to_vector(y_rep));
  int<lower=0, upper=1> p_mean = mean_y_rep >= mean_y;

  // posterior predictive sample std dev replicated successes
  real<lower=0> sd_y_rep = sd(to_vector(y_rep));
  int<lower=0, upper=1> p_sd = sd_y_rep >= sd_y;
}
