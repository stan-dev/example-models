data {
  int<lower=1> N;
  int<lower=1> J; // number of counties
  array[N] int<lower=1, upper=J> county;
  vector[N] u;
  vector[N] x;
  vector[N] y;
}
transformed data {
  matrix[N, 2] cov = [x', u']';
}
parameters {
  vector[2] beta;
  real mu_a;
  real<lower=0> sigma;
  real<lower=0> sigma_a;
  vector<offset=mu_a, multiplier=sigma_a>[J] alpha;
}
model {
  mu_a ~ std_normal();
  beta ~ normal(0, 100);
  sigma ~ cauchy(0, 2.5);
  sigma_a ~ cauchy(0, 2.5);
  alpha ~ normal(mu_a, sigma_a);
  
  y ~ normal_id_glm(cov, alpha[county], beta, sigma);
}
