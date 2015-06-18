## ---- irt-1pl-pin-stan ----
data {
  int<lower=0> I;
  int<lower=0> J;
  int<lower=0,upper=1> y[I,J];
}
parameters {
  vector[I] b;
  vector[J - 1] theta;
}
model {
  for (i in 1:I) {
    head(y[i], J - 1) ~ bernoulli_logit(theta - b[i]);
    y[i,J] ~ bernoulli_logit(b[i]);  // theta[J] = 0
  }
}
