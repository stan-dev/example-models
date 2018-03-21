/**
 * Bradley-Terry Model
 * Infer scores from incomplete set of paired comparisons.
 * Assumes ratings (a > b or b > a) are generated probabilistically
 * from scores.  Typical applications include those where observed
 * data are:
 *
 *  - winner of sporting match between two contestants (not all
 *  contestants need play each other and contestants might play each
 *  other more than once)
 *
 *  - whether product A is preferred to product B by a human rater
 *  (not all pairs of products need play each other and there may be
 *  multiple raters)
 *
 * Prior (i.e., population distribution of players [on log odds scale])
 *
 *   alpha[k] ~ normal(0, 3)
 *
 * Likelihood
 *
 *   y[n] ~ bernoulli(alpha[player1[n]] - alpha[player0[n]]))
 *
 * Pr[player 1 wins match n]
 *   = exp(alpha[player1[n]]) / (exp(alpha[player1[n]]) + exp(alpha[player0[n]]))
 *   = inverse_logit(alpha[player1[n]] - alpha[player0[n]])
 *
 * Odds for player 1 in match n
 *   = Pr[player 1 wins match n] / (1 - Pr[player 1 wins match n])
 *
 * Log odds for player 1 in match n
 *   = log(odds for player 1 in match n)
 *   = logit(Pr[player1 wins match n])
 *   = log(Pr[player 1 wins match n] / (1 - Pr[player 1 wins match n]))
 *
 * inverse_logit(u) = 1 / (1 + exp(-u)) = exp(u) / (1 + exp(u))
 *
 * Ranking
 *   player i better than player j = alpha[i] > alpha[j]
 *   best player = ARGMAX_k alpha[k]

 * Bradley, Ralph Allan; Terry, Milton E. 1952.  Rank analysis of
 *   incomplete block designs: I. The method of paired
 *   comparisons. Biometrika. 39 (3/4): 324. doi:10.2307/2334029
 */
data {
  int<lower = 0> K;                     // players
  int<lower = 0> N;                     // games
  int<lower=1, upper = K> player1[N];   // player 1 for game n
  int<lower=1, upper = K> player0[N];   // player 0 for game n
  int<lower = 0, upper = 1> y[N];       // winner for game n
}
parameters {
  vector[K] alpha;                      // ability for player n
}
model {
  alpha ~ normal(0, 1);
  y ~ bernoulli_logit(alpha[player1] - alpha[player0]);
}
generated quantities {
  int<lower=1, upper=K> ranking[K];       // rank of player ability
  {
    int ranked_index[K] = sort_indices_desc(alpha);
    for (k in 1:K)
      ranking[ranked_index[k]] = k;
  }
}
