//#  http://www.openbugs.net/Examples/Salm.html

//# the version without centering x's 
data {
  int<lower=0> Ndoses;
  int<lower=0> Nplates;
  array[Ndoses, Nplates] int<lower=0> y;
  array[Ndoses] real<lower=0> x;
}
parameters {
  real alpha;
  real beta;
  real gamma;
  real<lower=0> tau;
  array[Ndoses] vector[Nplates] lambda;
}
transformed parameters {
  real<lower=0> sigma;
  sigma = 1.0 / sqrt(tau);
}
model {
  alpha ~ normal(0.0, 100);
  beta ~ normal(0.0, 100);
  gamma ~ normal(0.0, 1.0E5);
  tau ~ gamma(0.001, 0.001);
  for (dose in 1 : Ndoses) {
    lambda[dose] ~ normal(0.0, sigma);
    y[dose] ~ poisson_log(alpha + beta * log(x[dose] + 10) + gamma * x[dose]
                          + lambda[dose]);
  }
}
