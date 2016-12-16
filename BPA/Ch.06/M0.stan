data {
  int<lower=0> M;               // Size of augumented data set
  int<lower=0> T;               // Number of sampling occasions
  int<lower=0,upper=1> y[M, T]; // Capture-history matrix
}

transformed data {
  int<lower=0> s[M];            // Totals in each row
  int<lower=0> C;               // Size of observed data set

  C = 0;
  for (i in 1:M) {
    s[i] = sum(y[i]);
    if (s[i] > 0)
      C = C + 1;
  }
}

parameters {
  real<lower=0,upper=1> omega;  // Inclusion probability
  real<lower=0,upper=1> p;      // Detection probability
}

model {
  // Priors are imlicitly defined;
  //  omega ~ uniform(0, 1);
  //  p ~ uniform(0, 1);

  // Likelihood
  for (i in 1:M) {
    real lp[2];

    if (s[i] > 0) {
      // z[i] == 1
      target += bernoulli_lpmf(1 | omega)
              + binomial_lpmf(s[i] | T, p);
    } else { // s[i] == 0
      // z[i] == 1
      lp[1] = bernoulli_lpmf(1 | omega)
            + binomial_lpmf(0 | T, p);
      // z[i] == 0
      lp[2] = bernoulli_lpmf(0 | omega);
      target += log_sum_exp(lp[1], lp[2]);
    }
  }
}

generated quantities {
  int<lower=C> N;

  N = C + binomial_rng(M, omega * (1 - p)^T);
}
