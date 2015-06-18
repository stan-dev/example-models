## ---- irt-1pl-stan ----
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
  theta ~ normal(0, 1);
  b ~ normal(-1, 2);
  for (i in 1:I)
    y[i] ~ bernoulli_logit(theta - b[i]);
}
