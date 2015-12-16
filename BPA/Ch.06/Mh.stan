data {
  int<lower=0> M;               // Size of augumented data set
  int<lower=0> T;               // Number of sampling occasions
  int<lower=0,upper=T> y[M];    // Capture-history matrix
}

transformed data {
  int<lower=0> C;               // Size of observed data set

  C <- 0;
  for (i in 1:M) {
    if (y[i] > 0)
      C <- C + 1;
  }
}

parameters {
  real<lower=0,upper=1> omega;  // Inclusion probability
  real<lower=0,upper=1> mean_p; // Detection probability
  real<lower=0> sigma;
  real<lower=-16,upper=16> eps[M];      // Individual random effects
}

model {
  real lp[2];

  // Priors
  omega ~ uniform(0, 1);        // Inclusion probability
  mean_p ~ uniform(0, 1);       // Mean detection probability
  sigma ~ uniform(0, 5);

  // Likelihood
  for (i in 1:M) {
    eps[i] ~ normal(logit(mean_p), sigma) T[-16, 16];
	     // See web appendix A in Royle (2009)

    if (y[i] > 0) {
      // z[i] == 1
      increment_log_prob(bernoulli_log(1, omega) +
                         binomial_logit_log(y[i], T, eps[i]));
    } else { // y[i] == 0
      // z[i] == 1
      lp[1] <- bernoulli_log(1, omega) +
               binomial_logit_log(0, T, eps[i]);
      // z[i] == 0
      lp[2] <- bernoulli_log(0, omega);
      increment_log_prob(log_sum_exp(lp));
    }
  }
}

generated quantities {
  int<lower=0,upper=1> zero[M];
  int<lower=C> N;

  for (i in 1:M) {
    real pr;

    pr <- pow(1.0 - inv_logit(eps[i]), T);
    zero[i] <- bernoulli_rng(omega * pr);
  }
  N <- C + sum(zero);
}
