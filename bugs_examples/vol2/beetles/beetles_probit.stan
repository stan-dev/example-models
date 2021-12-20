data {
  int<lower=0> N;
  array[N] int<lower=0> n;
  array[N] int<lower=0> r;
  vector[N] x;
}
transformed data {
  vector[N] centered_x;
  real mean_x;
  mean_x = mean(x);
  centered_x = x - mean_x;
}
parameters {
  real alpha_star;
  real beta;
}
transformed parameters {
  array[N] real p;
  for (i in 1 : N) {
    p[i] = Phi(alpha_star + beta * centered_x[i]);
  }
}
model {
  alpha_star ~ normal(0.0, 1.0);
  beta ~ normal(0.0, 1.0E4);
  r ~ binomial(n, p);
}
generated quantities {
  real alpha;
  array[N] real llike;
  array[N] real rhat;
  
  alpha = alpha_star - beta * mean_x;
  
  for (i in 1 : N) {
    llike[i] = r[i] * log(p[i]) + (n[i] - r[i]) * log(1 - p[i]);
    rhat[i] = p[i] * n[i];
  }
}
