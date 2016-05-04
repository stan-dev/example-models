// Binomial-mixture model with overdispersion in both abundance and detection
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
  vector[7] beta;
  vector[R] eps_raw;
  real<lower=0> sd_lam;
  real<lower=0> sd_p;
  real logit_p[R, T, 7];   // Originally `lp' in the BPA book
}

transformed parameters {
  vector[R] eps;
  matrix[R, 7] log_lambda;
  vector[K+1] lp[R, 7];

  eps <- sd_lam * eps_raw;
  for (i in 1:R) {                    // Loop over R sites (95)
    for (k in 1:7)
      log_lambda[i, k] <- alpha_lam[k] + eps[i];
    for (k in 1:(first[i] - 1))
      lp[i, k] <- rep_vector(negative_infinity(), K + 1);
    for (k in (last[i] + 1):7)
      lp[i, k] <- rep_vector(negative_infinity(), K + 1);
    for (k in first[i]:last[i]) {     // Loop over days
      lp[i, k, 1:max_y[i, k]] <- rep_vector(negative_infinity(), max_y[i, k]);
      for (n in max_y[i, k]:K)
        lp[i, k, n + 1] <- poisson_log_log(n, log_lambda[i, k])
          + binomial_logit_log(y[i, 1:T, k], n, logit_p[i, 1:T, k]);
    }
  }
}

model {
  // Priors
  alpha_lam ~ normal(0, sqrt(10));
  beta ~ normal(0, sqrt(10));
  eps_raw ~ normal(0, 1);
  // Half-Cauchy priors are used on sd_lam and sd_p, instead of
  // uniform(0, 3) used in the book;
  sd_lam ~ cauchy(0, 2.5);
  sd_p ~ cauchy(0, 2.5);

  for (i in 1:R)                      // Loop over R sites (95)
    for (k in first[i]:last[i])       // Loop over days
      logit_p[i, 1:T, k] ~ normal(beta[k], sd_p);

  // Likelihood
  for (i in 1:R)                      // Loop over R sites (95)
    for (k in first[i]:last[i])       // Loop over days
      increment_log_prob(log_sum_exp(lp[i, k, (max_y[i, k] + 1):(K + 1)]));
}

generated quantities {
  int totalN[7];
  vector[7] mean_abundance;
  vector[7] mean_N;
  vector[7] mean_detection;
  real fit;
  real fit_new;

  {
    int N[R, 7];
    real eval[R, T, 7];
    real y_new[R, T, 7];
    real p[R, T, 7];
    matrix[T, 7] E[R];
    matrix[T, 7] E_new[R];
    matrix[R, 7] ik_p;
    vector[7] num_obs_site;

    for (i in 1:R) {                         // Loop over R sites (95)
      for (k in 1:(first[i] - 1)) {
        N[i, k] <- 0;
        ik_p[i, k] <- 0;
        E[i, 1:T, k] <- rep_vector(0, T);
        E_new[i, 1:T, k] <- rep_vector(0, T);
      }
      for (k in (last[i] + 1):7) {
        N[i, k] <- 0;
        ik_p[i, k] <- 0;
        E[i, 1:T, k] <- rep_vector(0, T);
        E_new[i, 1:T, k] <- rep_vector(0, T);
      }
      for (k in first[i]:last[i]) {          // Loop over days (7)
        vector[K+1] pr;

        pr <- softmax(lp[i, k]);
        N[i, k] <- categorical_rng(pr) - 1;
        for (j in 1:T) {
          p[i, j, k] <- inv_logit(logit_p[i, j, k]);

          // Assess model fit using Chi-squared discrepancy
          // Compute fit statistic E for observed data
          eval[i, j, k] <- p[i, j, k] * N[i, k];   // Expected values
          E[i, j, k] <- square(y[i, j, k] - eval[i, j, k])
            / (eval[i, j, k] + 0.5);
          // Generate replicate data and compute fit stats for them
          y_new[i, j, k] <- binomial_rng(N[i, k], p[i, j, k]);
          E_new[i, j, k] <- square(y_new[i, j, k] - eval[i, j, k])
            / (eval[i, j, k] + 0.5);
        }
        ik_p[i, k] <- mean(p[i, 1:T, k]);
      }
    }

    for (k in 1:7) {
      num_obs_site[k] <- 0.0;
      for (i in 1:R)
        num_obs_site[k] <- num_obs_site[k] + (y[i, 1, k] != -1);

      totalN[k] <- sum(N[, k]);  // Total pop. size across all sites
      mean_abundance[k] <- mean(exp(log_lambda[, k]));
      mean_N[k] <- totalN[k] / num_obs_site[k];
      mean_detection[k] <- sum(ik_p[, k]) / num_obs_site[k];
     }
    fit <- 0.0;
    fit_new <- 0.0;
    for (i in 1:R) {
      fit <- fit + sum(E[i]);
      fit_new <- fit_new + sum(E_new[i]);
    }
  }
}
