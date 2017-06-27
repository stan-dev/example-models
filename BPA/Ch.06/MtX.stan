data {
  int<lower=0> M;               // Size of augumented data set
  int<lower=0> T;               // Number of sampling occasions
  int<lower=0> C;               // Size of observed data set
  int<lower=0,upper=1> y[M, T]; // Capture-history matrix
  real<lower=-6,upper=6> bsize[C];      // Body size
  real<lower=0> prior_sd_upper;
}

transformed data {
  int<lower=0> s[M];            // Totals in each row

  for (i in 1:M)
    s[i] = sum(y[i]);
}

parameters {
  real<lower=0,upper=1> omega;          // Inclusion probability
  vector<lower=0,upper=1>[T] mean_p;    // Mean detection probability
  real beta;
  real mu_size;
  real<lower=0,upper=prior_sd_upper> sd_size;
  real<lower=-6,upper=6> bsize_mis[M - C];      // Missing data
}

transformed parameters {
  vector[T] alpha = logit(mean_p);
  matrix[M, T] logit_p;

  for (i in 1:C)
    logit_p[i] = alpha' + beta * bsize[i];
  for (i in (C + 1):M)
    logit_p[i] = alpha' + beta * bsize_mis[i - C];
}

model {
  // Priors
  //  omega ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);
  beta ~ normal(0, 10);
  mu_size ~ normal(0, 10);
  //  sd_size ~ uniform(0, prior_sd_upper);   // Provide upper bound as data

  // Likelihood
  for (i in 1:C)
    bsize[i] ~ normal(mu_size, sd_size) T[-6, 6];
  for (i in (C + 1):M)
    bsize_mis[i - C] ~ normal(mu_size, sd_size) T[-6, 6];

  for (i in 1:M)
    if (s[i] > 0)
      // z[i] == 1
      target += bernoulli_lpmf(1 | omega)
              + bernoulli_logit_lpmf(y[i] | logit_p[i]);
    else // s[i] == 0
      target += log_sum_exp(bernoulli_lpmf(1 | omega)   // z[i] == 1
                            + bernoulli_logit_lpmf(0 | logit_p[i]),
                            bernoulli_lpmf(0 | omega)); // z[i] == 0
}

generated quantities {
  matrix<lower=0,upper=1>[M, T] p = inv_logit(logit_p);
  int<lower=0,upper=1> z[M];
  int<lower=C> N;

  for (i in 1:M) {
    if(s[i] > 0) {  // species was detected at least once
      z[i] = 1;
    } else {        // species was never detected
      // prob never detected given present
      real pr = prod(rep_vector(1, T) - p[i]');
      z[i] = bernoulli_rng(omega * pr / (omega * pr + (1 - omega)));
    }
  }
  N = sum(z);
}
