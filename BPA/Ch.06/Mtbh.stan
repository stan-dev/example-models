data {
  int<lower=0> M;               // Number of species
  int<lower=0> T;               // Number of grab occasions
  int<lower=0,upper=1> y[M, T]; // Grab counts
}

transformed data {
  int<lower=0> s[M];            // Detection times for each species
  int<lower=0> C;               // Number of observed species

  C = 0;
  for (i in 1:M) {
    s[i] = sum(y[i]);
    if (s[i] > 0)
      C = C + 1;
  }
}

parameters {
  real<lower=0,upper=1> omega;          // Inclusion probability
  vector<lower=0,upper=1>[T] mean_p;    // Mean detection probability
  real gamma;
  real<lower=0,upper=3> sigma;
  // In case a weakly informative prior is used
  //  real<lower=0> sigma;
  vector[M] eps_raw;
}

transformed parameters {
  vector[M] eps = sigma * eps_raw;
  vector[T] alpha = logit(mean_p);
  matrix[M, T] logit_p;

  // First occasion: no term for recapture (gamma)
  logit_p[ , 1] = alpha[1] + eps;

  // All subsequent occasions: includes recapture term (gamma)
  for (i in 1:M)
    for (j in 2:T)
      logit_p[i, j] = alpha[j] + eps[i] + gamma * y[i, j - 1];
}

model {
  // Priors
  gamma ~ normal(0, 10);
  // Uniform priors are implicitly defined.
  //  omega ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);
  //  sigma ~ uniform(0, 3);
  // In case a weakly informative prior is used
  //  sigma ~ normal(1.5, 0.75);
  eps_raw ~ normal(0, 1);

  // Likelihood
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
    if(s[i] > 0) {  // animal was detected at least once
      z[i] = 1;
    } else {        // animal never detected
      // prob never detected given present
      real pr = prod(rep_vector(1, T) - p[i]');
      z[i] = bernoulli_rng(omega * pr / (omega * pr + (1 - omega)));
    }
  }
  N = sum(z);
}
