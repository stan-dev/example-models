/* This program has a mistake in it, as will be explained later */

data {
  int N_teams;
  int N_games;
  vector[N_teams] prior_score;
  array[N_games] int team_1;
  array[N_games] int team_2;
  vector[N_games] score_1;
  vector[N_games] score_2;
  real df;
}
transformed data {
  vector[N_games] dif;
  vector[N_games] sqrt_dif;
  dif = score_1 - score_2;
  for (i in 1 : N_games) {
    sqrt_dif[i] = (step(dif[i]) - 0.5) * sqrt(abs(dif[i]));
  }
}
parameters {
  vector[N_teams] alpha;
  real b;
  real<lower=0> sigma_a;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[N_teams] a;
  a = b * prior_score + sigma_a * alpha;
}
model {
  alpha ~ normal(0, 1);
  sqrt_dif ~ student_t(df, a[team_1] - a[team_2], sigma_y);
}
