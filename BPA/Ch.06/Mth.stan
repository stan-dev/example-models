data {
  int<lower=0> M;               // Size of augumented data set
  int<lower=0> T;               // Number of sampling occasions
  int<lower=0,upper=1> y[M, T]; // Capture-history matrix
}

transformed data {
  int<lower=0,upper=T> s[M];    // Totals in each row
  int<lower=0,upper=M> C;       // Size of observed data set

  C = 0;
  for (i in 1:M) {
    s[i] = sum(y[i]);
    if (s[i] > 0)
      C = C + 1;
  }
}

parameters {
  real<lower=0,upper=1> omega;          // Inclusion probability
  real<lower=0,upper=1> mean_p[T];      // Mean detection probability
  real<lower=0,upper=5> sigma;
  // In case a weakly informative prior is used
  //  real<lower=0> sigma;
  vector[M] eps_raw;
}

transformed parameters {
  vector[M] eps = sigma * eps_raw;    // Random effects
  real mean_lp[T] = logit(mean_p);
  matrix[M, T] logit_p;

  for (j in 1:T)
    logit_p[ , j] = mean_lp[j] + eps;
}

model {
  // Priors are implicitly defined.
  //  omega ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);
  //  sigma ~ uniform(0, 5);
  // In case a weakly informative prior is used
  //  sigma ~ normal(2.5, 1.25);
  eps_raw ~ normal(0, 1);

  // Likelihood
  for (i in 1:M) {
    if (s[i] > 0)
      // z[i] == 1
      target += bernoulli_lpmf(1 | omega)
              + bernoulli_logit_lpmf(y[i] | logit_p[i]);
    else // s[i] == 0
      target += log_sum_exp(bernoulli_lpmf(1 | omega)   // z[i] == 1
                            + bernoulli_logit_lpmf(0 | logit_p[i]),
                            bernoulli_lpmf(0 | omega)); // z[i] == 0
  }
}

generated quantities {
  matrix<lower=0,upper=1>[M, T] p = inv_logit(logit_p);
  int<lower=0,upper=1> z[M];
  int<lower=C> N;

  for (i in 1:M) {
    if(s[i] > 0) {  // animal was detected at least once
      z[i] = 1;
    } else {        // animal was never detected
      // prob never detected given present
      real pr = prod(rep_vector(1, T) - p[i]');
      z[i] = bernoulli_rng(omega * pr / (omega * pr + (1 - omega)));
    }
  }
  N = sum(z);
}
