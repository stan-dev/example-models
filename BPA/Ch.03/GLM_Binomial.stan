data {
  int<lower=0> nyears;          // Number of Years
  int<lower=0> C[nyears];       // Counts
  int<lower=0> N[nyears];       // Binomial Totals
  real year[nyears];            // Year covariates
}

parameters {
  real alpha;
  real beta1;
  real beta2;
}

transformed parameters {
  real logit_p[nyears];

  for (i in 1:nyears)
    // link function and linear predictor
    logit_p[i] <- alpha +
                  beta1 * year[i] +
                  beta2 * pow(year[i], 2);
}

model {
  // Priors
  alpha ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  beta2 ~ normal(0, 100);

  // Likelihood
  for (i in 1:nyears) {
    // 1. Distribution for random part
    C[i] ~ binomial_logit(N[i], logit_p[i]);
  }
}

generated quantities {
  real<lower=0,upper=1> p[nyears];

  for (i in 1:nyears)
    p[i] <- inv_logit(logit_p[i]);
}
