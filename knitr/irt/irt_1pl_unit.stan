## ---- irt-1pl-unitt-stan ----
data {
  int<lower=0> I;
  int<lower=0> J;
  int<lower=0,upper=1> y[I,J];
}
parameters {
  vector[I] b;
  vector[J - 1] theta_raw;
}
transformed parameters {
  vector[J] theta;
  theta[1] <- -sum(theta_raw);
  for (j in 2:J)
    theta[j] <- theta_raw[j];
}
model {
  for (i in 1:I)
    y[i] ~ bernoulli_logit(theta - b[i]);
}
