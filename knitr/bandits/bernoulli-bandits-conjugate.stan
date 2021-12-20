/**
 *
 * Same model as bernoulli-bandits.stan with sufficient stats rather
 * than Bernoulli.
 */
data {
  int<lower=1> K; // num arms
  int<lower=0> N; // num trials
  array[N] int<lower=1, upper=K> z; // arm on trial n
  array[N] int<lower=0, upper=1> y; // reward on trial n
}
transformed data {
  array[K] int<lower=0> successes = rep_array(0, K);
  array[K] int<lower=0> trials = rep_array(0, K);
  for (n in 1 : N) {
    trials[z[n]] += 1;
    successes[z[n]] += y[n];
  }
}
generated quantities {
  simplex[K] is_best;
  vector<lower=0, upper=1>[K] theta;
  for (k in 1 : K) {
    theta[k] = beta_rng(1 + successes[k], 1 + trials[k] - successes[k]);
  }
  {
    real best_prob = max(theta);
    for (k in 1 : K) {
      is_best[k] = theta[k] >= best_prob;
    }
    is_best /= sum(is_best); // uniform for ties
  }
}
