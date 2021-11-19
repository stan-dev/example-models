data {
  int<lower=0> J; // number of schools
  array[J] real y; // estimated treatment effect (school j)
  array[J] real<lower=0> sigma; // std err of effect estimate (school j)
}
parameters {
  real mu;
  array[J] real theta;
  real<lower=0> tau;
}
model {
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
}
