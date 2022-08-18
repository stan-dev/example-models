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
  vector[J] alpha;
  vector[2] beta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] alphas = alpha[county];
}
model {
  y ~ normal_id_glm(xs, alphas, beta, sigma);  
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ normal(0, 10);
}
