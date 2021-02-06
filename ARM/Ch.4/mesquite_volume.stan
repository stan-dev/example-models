data {
  int<lower=0> N;
  vector[N] weight;
  vector[N] diam1;
  vector[N] diam2;
  vector[N] canopy_height;
}
transformed data {
  vector[N] log_weight = log(weight);
  vector[N] log_canopy_volume = log(diam1 .* diam2 .* canopy_height);
  matrix[N,1] x = [log_canopy_volume']';
}
parameters {
  real alpha;
  vector[1] beta;
  real<lower=0> sigma;
}
model {
  log_weight ~ normal_id_glm(x, alpha, beta, sigma);
}
