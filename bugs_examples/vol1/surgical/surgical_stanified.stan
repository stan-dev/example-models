# http://openbugs.net/Examples/Surgical.html
# random effects model
# Stanified: use sigma instead of squared, move p to generated quantitites
# Narrower priors on mu and sigma
data {
  int<lower=0> N;
  int r[N];
  int n[N];
}
parameters {
   real mu;
   real b[N];
   real<lower=0> sigma;
}

model {
  mu ~ normal(0.0, 20);
  sigma ~ cauchy(0,1); # remove sigma-squared
  b ~ normal(mu, sigma);
  r ~ binomial_logit(n, b);
}

generated quantities {
  real<lower=0,upper=1> p[N];
  real<lower=0> ranks[N];
  real pop_mean;
  pop_mean <- inv_logit(mu);
  for (i in 1:N){
    p[i] <- inv_logit(b[i]);
    ranks[i] <- rank(b,i);
  }
}
