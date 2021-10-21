data {
  real alpha;
  real beta;
  real<lower=0> sigma2;
  int<lower=0> J;
  array[J] int y;
  vector[J] Z;
  array[J] int n;
}
transformed data {
  real<lower=0> sigma;
  sigma = sqrt(sigma2);
}
parameters {
  real theta1;
  real theta2;
  vector[J] X;
}
model {
  array[J] real p;
  theta1 ~ normal(0, 32); // 32^2 = 1024 
  theta2 ~ normal(0, 32);
  X ~ normal(alpha + beta * Z, sigma);
  y ~ binomial_logit(n, theta1 + theta2 * X);
}
