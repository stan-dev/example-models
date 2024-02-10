/** include in generated quantities block of model with parameter
    vector<lower=0, upper=1>[N] theta; // chance of success 
    and data: N, K, y, K_new, y_new, min_y, max_y, mean_y, sd_y
*/

  // posterior predictive log density remaining trials
  real log_p_new = 0;
  for (n in 1 : N) {
    log_p_new = log_p_new + binomial_lpmf(y_new[n] | K_new[n], theta[n]);
  }

  // posterior predictions on new data 
  array[N] int<lower=0> z = binomial_rng(K_new, theta); 

  // Pr[some theta > 0.35]  
  int<lower=0, upper=1> some_ability_gt_350 = max(theta) > 0.35;
  // Pr[some player ability > 400]
  array[N] int<lower=0, upper=1> ability_gt_400;
  for (n in 1 : N) {
    ability_gt_400[n] = theta[n] > 0.400;
  }
  // Pr[some player season average > 400]
  array[N] int<lower=0, upper=1> avg_gt_400;
  for (n in 1 : N) {
    // y[n] - observed, z[n] - expected hits in rest of season
    avg_gt_400[n] = ((y[n] + z[n]) / (0.0 + K[n] + K_new[n])) > 0.400;
  }

  // replications for existing items  
  array[N] int<lower=0> y_rep = binomial_rng(K, theta);
  
  // posterior predictive min replicated successes, test stat p_val
  real<lower=0> min_y_rep = min(y_rep);
  int<lower=0, upper=1> p_min = min_y_rep >= min_y;

  // posterior predictive max replicated successes, test stat p_val
  real<lower=0> max_y_rep = max(y_rep); 
  int<lower=0, upper=1> p_max = max_y_rep >= max_y;

  // posterior predictive sample mean replicated successes, test stat p_val
  real<lower=0> mean_y_rep = mean(to_vector(y_rep));
  int<lower=0, upper=1> p_mean = mean_y_rep >= mean_y;

  // posterior predictive sample std dev replicated successes, test stat p_val
  real<lower=0> sd_y_rep = sd(to_vector(y_rep));
  int<lower=0, upper=1> p_sd = sd_y_rep >= sd_y;
