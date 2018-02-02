/**
 * Straightforward Bernoulli formulation.
 */
data {
  int<lower = 1> K;                       // num arms
  int<lower = 0> N;                       // num trials
  int<lower = 1, upper = K> z[N];         // arm on trial n
  int<lower = 0, upper = 1> y[N];         // reward on trial n
}
parameters {
  vector<lower = 0, upper = 1>[K] theta;  // arm return prob
}
model {
  y ~ bernoulli(theta[z]);                // i.i.d. by arm
}
generated quantities {
  simplex[K] is_best;  // one hot or uniform with ties
  {
    real best_prob = max(theta);
    for (k in 1:K)
      is_best[k] = (theta[k] >= best_prob);
    is_best /= sum(is_best);  // uniform for ties
  }
}
