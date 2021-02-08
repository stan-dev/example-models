data {
  int<lower=1> N;
  int<lower=1> J; # number of counties
  int<lower=1,upper=J> county[N];
  vector[N] u;
  vector[N] x;
  vector[N] y;
}
transformed data {
  matrix[N,2] cov = [x', u']';
}
parameters {
  real mu_alpha;
  real mu_beta;
  real<lower=0> sigma;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  vector<offset=mu_alpha, multiplier=sigma_alpha>[J] alpha;
  vector<offset=mu_beta, multiplier=sigma_beta>[2] beta;
}
model {
  alpha ~ normal(mu_alpha, sigma_alpha);
  beta ~ normal(mu_beta, sigma_beta);
  sigma ~ cauchy(0, 2.5);
  mu_alpha ~ std_normal();
  sigma_alpha ~ cauchy(0, 2.5);
  mu_beta ~ std_normal();
  sigma_beta ~ cauchy(0, 2.5);

  y ~ normal_id_glm(cov, alpha[county], beta, sigma);
}
