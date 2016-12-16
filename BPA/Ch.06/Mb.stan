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
  real<lower=0,upper=1> p;      // Capture probability
                                //(not captureed during the preceeding occasion)
  real<lower=0,upper=1> c;      // Capture probability
                                //(captureed during the preceeding occasion)
}

transformed parameters {
  vector<lower=0,upper=1>[T] p_eff[M];

  for (i in 1:M) {
    // First occasion
    p_eff[i][1] = p;

    // All subsequent occasions
    for (j in 2:T)
      p_eff[i][j] = (1 - y[i, j - 1]) * p + y[i, j - 1] * c;
  }
}

model {
  // Priors are implicitly define;
  //  omega ~ uniform(0, 1);
  //  p ~ uniform(0, 1);

  // Likelihood
  for (i in 1:M) {
    real lp[2];

    if (s[i] > 0) {
      // z[i] == 1
      target += bernoulli_lpmf(1 |  omega)
              + bernoulli_lpmf(y[i] |  p_eff[i]);
    } else { // s[i] == 0
      // z[i] == 1
      lp[1] = bernoulli_lpmf(1 |  omega)
            + bernoulli_lpmf(0 | p_eff[i]);
      // z[i] == 0
      lp[2] = bernoulli_lpmf(0 | omega);
      target += log_sum_exp(lp[1], lp[2]);
    }
  }
}

generated quantities {
  int<lower=0,upper=1> zero[M];
  int<lower=C> N;
  real trap_response;

  for (i in 1:M) {
    real pr;

    pr = prod(rep_vector(1.0, T) - p_eff[i]);
    zero[i] = bernoulli_rng(omega * pr);
  }
  N = C + sum(zero);
  trap_response = c - p;
}
