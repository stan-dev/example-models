data {
  int<lower=0> M;
  int<lower=0> T;
  int<lower=0,upper=1> y[M, T];
}

transformed data {
  int<lower=0,upper=T> s[M];
  int<lower=0,upper=M> C;

  C <- 0;
  for (i in 1:M) {
    s[i] <- sum(y[i]);
    if (s[i] > 0)
      C <- C + 1;
  }
}

parameters {
  real<lower=0,upper=1> omega;
  real<lower=0,upper=1> mean_p[T];
  real<lower=0> sigma;
  real eps[M];
}

transformed parameters {
  real mean_lp[T];
  vector<lower=0,upper=1>[T] p[M];

  for (j in 1:T)
    mean_lp[j] <- logit(mean_p[j]); // Define logit
  for (i in 1:M)
    for (j in 1:T)
      p[i][j] <- inv_logit(mean_lp[j] + eps[i]);
}

model {
  real lp[2];

  // Priors
  omega ~ uniform(0, 1);
  mean_p ~ uniform(0, 1);
  sigma ~ uniform(0, 5);

  // Likelihood
  for (i in 1:M) {
    eps[i] ~ normal(0, sigma) T[-16, 16];
             // See web appendix A in Royle (2009)

    if (s[i] > 0) {
      // z[i] == 1
      increment_log_prob(bernoulli_log(1, omega) +
                         bernoulli_log(y[i], p[i]));
    } else { // s[i] == 0
      // z[i] == 1
      lp[1] <- bernoulli_log(1, omega) +
               bernoulli_log(0, p[i]);
      // z[i] == 0
      lp[2] <- bernoulli_log(0, omega);
      increment_log_prob(log_sum_exp(lp));
    }
  }
}

generated quantities {
  int<lower=0,upper=1> zero[M];
  int<lower=C> N;

  for (i in 1:M) {
    real pr;

    pr <- prod(rep_vector(1.0, T) - p[i]);
    zero[i] <- bernoulli_rng(omega * pr);
  }
  N <- C + sum(zero);
}
