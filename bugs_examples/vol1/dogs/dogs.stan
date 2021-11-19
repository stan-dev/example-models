// compared with JAGS version in 
// the R package BUGSExamples (https://r-forge.r-project.org/R/?group_id=882) 
data {
  int<lower=0> Ndogs;
  int<lower=0> Ntrials;
  array[Ndogs, Ntrials] int Y;
}
transformed data {
  array[Ndogs, Ntrials] int y;
  array[Ndogs, Ntrials] int xa;
  array[Ndogs, Ntrials] int xs;
  for (dog in 1 : Ndogs) {
    xa[dog, 1] = 0;
    xs[dog, 1] = 0;
    for (trial in 2 : Ntrials) {
      for (k in 1 : (trial - 1)) {
        xa[dog, trial] = xa[dog, trial] + Y[dog, k];
      }
      xs[dog, trial] = trial - 1 - xa[dog, trial];
    }
  }
  for (dog in 1 : Ndogs) {
    for (trial in 1 : Ntrials) {
      y[dog, trial] = 1 - Y[dog, trial];
    }
  }
}
parameters {
  real<upper=-0.00001> alpha;
  real<upper=-0.00001> beta;
}
model {
  alpha ~ normal(0.0, 316.2);
  beta ~ normal(0.0, 316.2);
  for (dog in 1 : Ndogs) {
    for (trial in 2 : Ntrials) {
      y[dog, trial] ~ bernoulli(exp(alpha * xa[dog, trial]
                                    + beta * xs[dog, trial]));
    }
  }
}
generated quantities {
  real A;
  real B;
  A = exp(alpha);
  B = exp(beta);
}
