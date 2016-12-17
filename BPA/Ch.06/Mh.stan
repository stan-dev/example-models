data {
  int<lower=0> M;               // Size of augumented data set
  int<lower=0> T;               // Number of sampling occasions
  int<lower=0,upper=T> y[M];    // Capture-history matrix
}

transformed data {
  int<lower=0> C;               // Size of observed data set

  C = 0;
  for (i in 1:M) {
    if (y[i] > 0)
      C = C + 1;
  }
}

parameters {
  real<lower=0,upper=1> omega;      // Inclusion probability
  real<lower=0,upper=1> mean_p;     // Mean detection probability
  real<lower=0> sigma;
  real<lower=-16,upper=16> eps[M];  // Individual random effects
}

model {
  real lp[2];

  // Priors are implicitly defined.
  //  omega ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);
  // An improper uniform prior is used for sigma
  //  instead of uniform(0, 5).

  // Likelihood
  for (i in 1:M) {
    eps[i] ~ normal(logit(mean_p), sigma) T[-16, 16];
	     // See web appendix A in Royle (2009)

    if (y[i] > 0) {
      // z[i] == 1
      target += bernoulli_lpmf(1 | omega)
              + binomial_logit_lpmf(y[i] | T, eps[i]);
    } else { // y[i] == 0
      // z[i] == 1
      lp[1] = bernoulli_lpmf(1 | omega)
            + binomial_logit_lpmf(0 | T, eps[i]);
      // z[i] == 0
      lp[2] = bernoulli_lpmf(0 | omega);
      target += log_sum_exp(lp[1], lp[2]);
    }
  }
}

generated quantities {
  int<lower=0,upper=1> zero[M];
  int<lower=C> N;

  for (i in 1:M) {
    real pr;

    pr = (1.0 - inv_logit(eps[i]))^T;
    zero[i] = bernoulli_rng(omega * pr);
  }
  N = C + sum(zero);
}
