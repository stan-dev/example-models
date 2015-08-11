data {
  int<lower=0> I;
  int<lower=0> J;
  int<lower=0,upper=1> y[I,J];
}
parameters {
  vector[I] b;
  vector[J] theta;
}
model {
  for (i in 1:I)
    y[i] ~ bernoulli_logit(theta - b[i]);
}
