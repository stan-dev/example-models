data {
  int<lower=1> N;
  int<lower=1> J; // number of counties
  array[N] int<lower=1, upper=J> county;
  vector[N] x;
  vector[N] y;
}
transformed data {
  matrix[N, 1] cov = [x']';
}
parameters {
  vector[1] beta;
  real<lower=0> sigma_a;
  real<lower=0> sigma_y;
  real mu_a;
  vector<offset=mu_a, multiplier=sigma_a>[J] a;
}
model {
  beta ~ std_normal();
  mu_a ~ std_normal();
  sigma_a ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  
  a ~ normal(mu_a, sigma_a);
  y ~ normal_id_glm(cov, a[county], beta, sigma_y);
}
