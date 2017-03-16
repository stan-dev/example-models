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
  real<lower=0,upper=1> mean_p[T];      // Mean detection probability
  real beta;
  real mu_size;
  real<lower=0,upper=prior_sd_upper> sd_size;
  real<lower=-6,upper=6> bsize_mis[M-C];        // Missing data
}

transformed parameters {
  real alpha[T];
  vector<lower=0,upper=1>[T] p[M];

  for (j in 1:T)
    alpha[j] = logit(mean_p[j]); // Define logit
  for (i in 1:C)
    for (j in 1:T)
      p[i][j] = inv_logit(alpha[j] + beta * bsize[i]);
  for (i in (C + 1):M)
    for (j in 1:T)
      p[i][j] = inv_logit(alpha[j] + beta * bsize_mis[i - C]);
}

model {
  real lp[2];

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

  for (i in 1:M) {  // Loop over individuals
    if (s[i] > 0) {
      // z[i] == 1
      target += bernoulli_lpmf(1 | omega)
              + bernoulli_lpmf(y[i] | p[i]);
    } else { // s[i] == 0
      // z[i] == 1
      lp[1] = bernoulli_lpmf(1 | omega)
            + bernoulli_lpmf(0 | p[i]);
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

    pr = prod(rep_vector(1.0, T) - p[i]);
    zero[i] = bernoulli_rng(omega * pr);
  }
  N = C + sum(zero);
}
