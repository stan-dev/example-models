/** observed outcome y for N items in K, K_new binary trials **/
data {
  int<lower=0> N; // items
  array[N] int<lower=0> K; // initial trials
  array[N] int<lower=0> y; // initial successes
  
  array[N] int<lower=0> K_new; // new trials
  array[N] int<lower=0> y_new; // new successes
}
transformed data {
  // summary stats for observed data
  real min_y = min(y);
  real max_y = max(y);
  real mean_y = mean(to_vector(y));
  real sd_y = sd(to_vector(y));
}
