data {
  int<lower=0> N;                     // players
  int<lower=0> K;                     // first at-bats
  int<lower=0, upper=K> y[N];         // first hits
  int<lower=0> J[N];                  // remaining at bats
}
parameters {
  real<lower=0, upper=1> phi;         // population ability
  real<lower=1> kappa;                // population concentration
  vector<lower=0, upper=1>[N] theta;  // player ability
}
model {
  kappa ~ pareto(1, 1.5);                        // hyperprior
  theta ~ beta(phi * kappa, (1 - phi) * kappa);  // prior
  y ~ binomial(K, theta);                        // likelihood
}
generated quantities {
  int<lower=0, upper=1> avg_gt_400[N];      // season avg >= 0.400
  int<lower=0, upper=1> ability_gt_400[N];  // ability >= 0.400
  for (n in 1:N) {
    int z_n;
    z_n <- binomial_rng(J[n], theta[n]);
    avg_gt_400[n] <- (((y[n] + z_n) / (0.0 + K + J[n])) > 0.400);
    ability_gt_400[n] <- (theta[n] > 0.400);
  }
}
