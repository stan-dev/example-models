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
  real<lower=0,upper=1> mean_p[T];      // Mean detection probability
  real gamma;
  real<lower=0,upper=3> sigma;
  // In case a weakly informative prior is used
  //  real<lower=0> sigma;
  real<lower=-16,upper=16> eps[M];
}

transformed parameters {
  real alpha[T];
  vector<lower=0,upper=1>[T] p[M];

  for (j in 1:T)
    alpha[j] = logit(mean_p[j]); // Define logit
  for (i in 1:M) {
    real logit_p;

    // First occasion: no term for recapture (gamma)
    logit_p = alpha[1] + eps[i];
    p[i][1] = inv_logit(logit_p);

    // All subsequent occasions: includes recapture term (gamma)
    for (j in 2:T) {
      logit_p = alpha[j] + eps[i] + gamma * y[i, j - 1];
      p[i][j] = inv_logit(logit_p);
    }
  }
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

  // Likelihood
  for (i in 1:M) {
    real lp[2];

    eps[i] ~ normal(0, sigma) T[-16, 16];
             // See web appendix A in Royle (2009)

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

  for (i in 1:M) {
    real pr;


    pr = prod(rep_vector(1.0, T) - p[i]);
    zero[i] = bernoulli_rng(omega * pr);
  }
  N = C + sum(zero);
}
