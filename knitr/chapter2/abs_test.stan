data {
  int N;
  vector[N] y;
}
transformed data {
  vector[N] abs_y = fabs(y);
}
parameters {
  real theta;
}
model {
  y ~ normal(theta, 1);
}
