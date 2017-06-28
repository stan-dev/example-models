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
  for (i in 1:M)
    if (s[i] > 0)
      // z[i] == 1
      target += bernoulli_lpmf(1 | omega)
              + binomial_lpmf(s[i] | T, p);
    else // s[i] == 0
      target += log_sum_exp( bernoulli_lpmf(1 | omega)   // z[i] == 1
                             + binomial_lpmf(0 | T, p),
                             bernoulli_lpmf(0 | omega)); // z[i] == 0
}

generated quantities {
  // prob present given never detected
  real omega_nd = (omega * (1 - p)^T) / (omega * (1 - p)^T + (1 - omega));
  int<lower=C, upper=M> N = C + binomial_rng(M - C, omega_nd);
}
