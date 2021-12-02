data {
  int<lower=0> I;
  int<lower=0> J;
  array[I, J] int<lower=0, upper=1> y;
}
parameters {
  vector[I] b;
  vector[J] theta;
}
model {
  for (i in 1 : I) {
    y[i] ~ bernoulli_logit(theta - b[i]);
  }
}
