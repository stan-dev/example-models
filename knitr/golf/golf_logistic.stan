data {
  int J;
  array[J] int n;
  vector[J] x;
  array[J] int y;
}
parameters {
  real a;
  real b;
}
model {
  y ~ binomial_logit(n, a + b * x);
}
