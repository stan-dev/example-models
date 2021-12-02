data {
  int<lower=1> N;
  int<lower=1> J; // number of counties
  array[N] int<lower=1, upper=J> county;
  vector[N] y;
}
parameters {
  real mu_a;
  real<lower=0> sigma_a;
  real<lower=0> sigma_y;
  vector<offset=(10 * mu_a), multiplier=sigma_a>[J] alpha;
}
model {
  mu_a ~ std_normal();
  alpha ~ normal(10 * mu_a, sigma_a);
  sigma_a ~ cauchy(0, 2.5);
  sigma_y ~ cauchy(0, 2.5);
  
  y ~ normal(alpha[county], sigma_y);
}
