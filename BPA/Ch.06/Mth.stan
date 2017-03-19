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
  real eps[M];                          // Random effects
}

transformed parameters {
  real mean_lp[T];
  vector<lower=0,upper=1>[T] p[M];

  for (j in 1:T)
    mean_lp[j] = logit(mean_p[j]); // Define logit
  for (i in 1:M)
    for (j in 1:T)
      p[i][j] = inv_logit(mean_lp[j] + eps[i]);
}

model {
  // Priors are implicitly defined.
  //  omega ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);
  //  sigma ~ uniform(0, 5);
  // In case a weakly informative prior is used
  //  sigma ~ normal(2.5, 1.25);

  // Likelihood
  for (i in 1:M) {
    eps[i] ~ normal(0, sigma) T[-16, 16];
             // See web appendix A in Royle (2009)

    if (s[i] > 0) {
      // z[i] == 1
      target += bernoulli_lpmf(1 | omega)
              + bernoulli_lpmf(y[i] | p[i]);
    } else { // s[i] == 0
      real lp[2];

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
  int<lower=0,upper=1> z[M];
  int<lower=C> N;

  for (i in 1:M) {
    if(s[i] > 0) {  // animal was detected at least once
      z[i] = 1;
    } else {        // animal was never detected
      real pr;      // prob never detected given present
      pr = prod(rep_vector(1.0, T) - p[i]);
      z[i] = bernoulli_rng(omega * pr / (omega * pr + (1 - omega)));
    }
  }
  N = sum(z);}
