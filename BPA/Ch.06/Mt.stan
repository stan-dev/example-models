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
  vector<lower=0,upper=1>[T] p; // Detection probability
}

model {
  // Priors are implicitly defined.
  //  omega ~ uniform(0, 1);
  //  p ~ uniform(0, 1);

  // Likelihood
  for (i in 1:M) {
    real lp[2];

    if (s[i] > 0) {
      // z[i] == 1
      target += bernoulli_lpmf(1 | omega)
              + bernoulli_lpmf(y[i] | p);
    } else { // s[i] == 0
      // z[i] == 1
      lp[1] = bernoulli_lpmf(1 | omega)
            + bernoulli_lpmf(0 | p);
      // z[i] == 0
      lp[2] = bernoulli_lpmf(0 | omega);
      target += log_sum_exp(lp[1], lp[2]);
    }
  }
}

generated quantities {
  int<lower=C> N;
  real pr;        // prob never captured given present
  real omega_nd;  // prob present given never captured; same for all animals
  
  pr = prod(rep_vector(1.0, T) - p);
  omega_nd = (omega * pr) / (omega * pr + (1 - omega));
  N = C + binomial_rng(M - C, omega_nd);
}
