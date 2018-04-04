/**
 * Bradley-Terry model for maximum likelihood estimation (i.e., no prior).
 */
data {
  int<lower = 0> K;                     // players
  int<lower = 0> N;                     // games
  int<lower=1, upper = K> player1[N];   // player 1 for game n
  int<lower=1, upper = K> player0[N];   // player 0 for game n
  int<lower = 0, upper = 1> y[N];       // winner for game n
}
parameters {
  vector[K - 1] alpha_raw;              // ability for players 1:K-1
}
transformed parameters {
  // enforces sum(alpha) = 0 for identifiability
  vector[K] alpha = append_row(alpha_raw, -sum(alpha_raw));
}
model {
  y ~ bernoulli_logit(alpha[player1] - alpha[player0]);
}
generated quantities {
  int<lower=1, upper=K> ranked[K] = sort_indices_desc(alpha);
}
