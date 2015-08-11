## ---- irt-1pl-adjust-stan ----
data {
  int<lower=0> I;
  int<lower=0> J;
  int<lower=0,upper=1> y[I,J];
}
parameters {
  vector[I] b_raw;
  vector[J] theta_raw;
}
transformed parameters {
  vector[I] b;
  vector[J] theta;
  { 
    real mean_theta_raw;
    mean_theta_raw <- mean(theta_raw);
    theta <- theta_raw - mean_theta_raw;
    b <- b_raw - mean_theta_raw;
  }
}
model {
  for (i in 1:I)
    y[i] ~ bernoulli_logit(theta - b[i]);
}
