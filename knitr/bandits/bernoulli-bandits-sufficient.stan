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
  array[2] int<lower=0> successes = rep_array(0, K);
  array[2] int<lower=0> trials = rep_array(0, K);
  for (n in 1 : N) {
    trials[z[n]] += 1;
    successes[z[n]] += y[n];
  }
}
parameters {
  vector<lower=0, upper=1>[K] theta; // arm return prob
}
model {
  successes ~ binomial(trials, theta);
}
generated quantities {
  simplex[K] is_best; // uniform over max
  {
    real best_prob = max(theta);
    for (k in 1 : K) {
      is_best[k] = theta[k] >= best_prob;
    }
    is_best /= sum(is_best); // uniform for ties
  }
}
