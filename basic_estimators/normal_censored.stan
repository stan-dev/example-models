data {
  real U;
  int<lower=0> N_censored;
  int<lower=0> N_observed;
  real<upper=U> y[N_observed];
}
parameters { 
  real mu;
}
model {
  for (n in 1:N_observed)
    y[n] ~ normal(mu,1.0) T[,U];
  target += N_censored * normal_lccdf(U | mu, 1);
}


