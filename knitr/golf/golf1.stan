data {
  int J;
  array[J] int n;
  vector[J] x;
  array[J] int y;
  real r;
  real R;
}
transformed data {
  vector[J] threshold_angle = asin((R - r) ./ x);
}
parameters {
  real<lower=0> sigma;
}
model {
  vector[J] p = 2 * Phi(threshold_angle / sigma) - 1;
  y ~ binomial(n, p);
}
generated quantities {
  real sigma_degrees = (180 / pi()) * sigma;
}
