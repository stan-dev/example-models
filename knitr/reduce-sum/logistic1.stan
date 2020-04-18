functions {
  real partial_sum(int[] slice_n_redcards,
                   int start, int end,
                   int[] n_games,
                   vector rating,
                   vector beta) {
    return binomial_logit_lpmf(slice_n_redcards |
                               n_games[start:end],
                               beta[1] + beta[2] * rating[start:end]);
  }
}
data {
  int N;
  int n_redcards[N];
  int n_games[N];
  vector[N] rating;
}
parameters {
  vector[2] beta;
}
model {
  int grainsize = 1;

  beta[1] ~ normal(0, 10);
  beta[2] ~ normal(0, 1);

  target += reduce_sum(partial_sum, n_redcards, grainsize,
                       n_games, rating, beta);
}