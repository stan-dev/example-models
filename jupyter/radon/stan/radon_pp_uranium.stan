data {
  int<lower=1> N;  // observations
  int<lower=1> J;  // counties
  array[N] int<lower=1, upper=J> county;
  vector[N] y;     // radon
  vector[N] x;     // floor
  vector[N] u;     // uranium
}
transformed data {
  matrix[N, 2] xs = [x', u']';
}
parameters {
  real mu_alpha;
  real<lower=0> sigma_alpha;
  vector<offset=mu_alpha, multiplier=sigma_alpha>[J] alpha;  // non-centered parameterization
  vector[2] beta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] alphas = alpha[county];
}
model {
  y ~ normal_id_glm(xs, alphas, beta, sigma);  
  alpha ~ normal(mu_alpha, sigma_alpha); // partial-pooling
  beta ~ normal(0, 10);
  sigma ~ normal(0, 10);
  mu_alpha ~ normal(0, 10);
  sigma_alpha ~ normal(0, 10);
}
generated quantities {
  array[N] real y_rep = normal_rng(alpha[county] + beta * xs, sigma);
}
