// Open-population binomial-mixuture model
data {
  int<lower=1> R;                // Number of sites
  int<lower=1> T;                // Number of replications
  int<lower=-1> y[R, T, 7];      // Counts (-1:NA)
  int<lower=1,upper=7> first[R]; // First occasion
  int<lower=1,upper=7> last[R];  // Last occasion
  int<lower=0> K;                // Upper bounds of population size
}

transformed data {
  int<lower=0> max_y[R, 7];

  for (i in 1:R) {
    for (k in 1:(first[i] - 1))
      max_y[i, k] <- 0;
    for (k in (last[i] + 1):7)
      max_y[i, k] <- 0;
    for (k in first[i]:last[i])
      max_y[i, k] <- max(y[i, 1:T, k]);
  }
}

parameters {
  vector[7] alpha_lam;
  vector<lower=0,upper=1>[7] p;  // Capture probability
}

transformed parameters {
  vector[K+1] lp[R, 7];          // Log probability

  // Likelihood
  for (i in 1:R) {               // Loop over R sites (95)
    for (k in 1:(first[i] - 1)) {
      lp[i, k, 1:(K + 1)] <- rep_vector(negative_infinity(), K + 1);
    }
    for (k in (last[i] + 1):7) {
      lp[i, k, 1:(K + 1)] <- rep_vector(negative_infinity(), K + 1);
    }
    for (k in first[i]:last[i]) {         // Loop over days
      lp[i, k, 1:max_y[i, k]] <- rep_vector(negative_infinity(), max_y[i, k]);
      for (n in max_y[i, k]:K) {
        lp[i, k, n + 1] <- poisson_log_log(n, alpha_lam[k])
          + binomial_log(y[i, 1:T, k], n, p[k]);
      }
    }
  }
}

model {
  // Priors
  alpha_lam ~ normal(0, 10);
  // A flat prior Uniform(0, 1) is implicitly used on p.

  // Likelihood
  for (i in 1:R)
    for (k in first[i]:last[i])
      increment_log_prob(log_sum_exp(lp[i, k, (max_y[i, k] + 1):(K + 1)]));
}

generated quantities {
  int totalN[7];
  real fit;
  real fit_new;
  vector[7] mean_abundance;

  {
    int N[R, 7];
    real eval[R, 7];
    real y_new[R, T, 7];
    matrix[T, 7] E[R];
    matrix[T, 7] E_new[R];

    for (i in 1:R) {                         // Loop over R sites (95)
      // Initialize N, E and E_new
      for (k in 1:(first[i] - 1))
        N[i, k] <- 0;
      for (k in (last[i] + 1):7)
        N[i, k] <- 0;
      E[i] <- rep_matrix(0, T, 7);
      E_new[i] <- rep_matrix(0, T, 7);

      for (k in first[i]:last[i]) {          // Loop over days
        vector[K+1] pr;

        pr <- softmax(lp[i, k]);
        N[i, k] <- categorical_rng(pr) - 1;
        eval[i, k] <- p[k] * N[i, k];        // Expected values
        for (j in 1:T) {
          // Assess model fit using Chi-squared discrepancy
          // Compute fit statistic E for observed data
          E[i, j, k] <- square(y[i, j, k] - eval[i, k])
            / (eval[i, k] + 0.5);
          // Generate replicate data and
          // Compute fit statistic E_new for replicate data
          y_new[i, j, k] <- binomial_rng(N[i, k], p[k]);
          E_new[i, j, k] <- square(y_new[i, j, k] - eval[i, k])
            / (eval[i, k] + 0.5);
        }
      }
    } //i
    for (k in 1:7)
      totalN[k] <- sum(N[1:R, k]);  // Total pop. size across all sites
    fit <- 0.0;
    fit_new <- 0.0;
    for (i in 1:R) {
      fit <- fit + sum(E[i]);
      fit_new <- fit_new + sum(E_new[i]);
    }
  }
  mean_abundance <- exp(alpha_lam);
}
