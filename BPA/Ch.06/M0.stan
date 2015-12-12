data {
  int<lower=0> M;
  int<lower=0> T;
  int<lower=0,upper=1> y[M, T];
}

transformed data {
  int<lower=0> s[M];
  int<lower=0> C;

  C <- 0;
  for (i in 1:M) {
    s[i] <- sum(y[i]);
    if (s[i] > 0)
      C <- C + 1;
  }
}

parameters {
  real<lower=0,upper=1> omega;
  real<lower=0,upper=1> p;
}

model {
  // Priors
  omega ~ uniform(0, 1);
  p ~ uniform(0, 1);

  // Likelihood
  for (i in 1:M) {
    real lp[2];

    if (s[i] > 0) {
      // z[i] == 1
      increment_log_prob(bernoulli_log(1, omega) +
                         binomial_log(s[i], T, p));
    } else { // s[i] == 0
      // z[i] == 1
      lp[1] <- bernoulli_log(1, omega) +
               binomial_log(0, T, p);
      // z[i] == 0
      lp[2] <- bernoulli_log(0, omega);
      increment_log_prob(log_sum_exp(lp));
    }
  }
}

generated quantities {
  int<lower=C> N;

  N <- C + binomial_rng(M, omega * pow(1 - p, T));
}
