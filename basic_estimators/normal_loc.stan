transformed data {
  vector[5] y = [2, 1, -0.5, 3, 0.25]';
}
parameters {
  real mu;
}
model {
  y ~ normal(mu, 1);
}
