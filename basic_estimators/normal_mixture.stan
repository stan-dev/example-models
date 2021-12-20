// normal mixture, unknown proportion and means, known variance
// p(y|mu,theta) = theta * Normal(y|mu[1],1) + (1-theta) * Normal(y|mu[2],1);

data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real<lower=0, upper=1> theta;
  vector[2] mu;
}
model {
  theta ~ uniform(0, 1); // equivalently, ~ beta(1,1);
  mu ~ normal(0, 10);
  for (n in 1 : N) {
    target += log_mix(theta, normal_lupdf(y[n] | mu[1], 1),
                      normal_lupdf(y[n] | mu[2], 1));
  }
}
