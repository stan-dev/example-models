// Binomial-mixture model with overdispersion in both abundance and detection

functions {
  /**
   * Returns log likelihood of N-mixture model
   * with 2 replicated observations using
   * bivariate Poisson distibution
   *
   * References
   * Dennis et al. (2015) Computational aspects of N-mixture models.
   *   Biometrics 71:237--246. DOI:10.1111/biom.12246
   * Stan users mailing list
   *   https://groups.google.com/forum/#!topic/stan-users/9mMsp1oB69g
   *
   * @param n          Number of observed individuals
   * @param log_lambda Log of Poisson mean of population size
   * @param logit_p    Logit of detection probability
   *
   * return Log probability
  */
  real bivariate_poisson_log_lpmf(int[] n, real log_lambda, real[] p) {
    real s[min(n) + 1];

    if (size(n) != 2)
      reject("Size of n must be 2.");
    if (p[1] < 0 || p[1] > 1 || p[2] < 0 || p[2] > 1)
      reject("p must be in [0,1].");
    for (u in 0:min(n))
      s[u + 1] = poisson_log_lpmf(n[1] - u | log_lambda + log(p[1]) + log1m(p[2]))
        + poisson_log_lpmf(n[2] - u | log_lambda + log(p[2]) + log1m(p[1]))
        + poisson_log_lpmf(u | log_lambda + log(p[1]) + log(p[2]));
    return log_sum_exp(s);
  }
}

data {
  int<lower=1> R;                // Number of sites
  int<lower=1> T;                // Number of replications; fixed as 2
  int<lower=-1> y[R, 2, 7];      // Counts (-1:NA)
  int<lower=1,upper=7> first[R]; // First occasion
  int<lower=1,upper=7> last[R];  // Last occasion
  int<lower=0> K;                // Upper bounds of population size
}

transformed data {
  int<lower=0> max_y[R, 7];
  int<lower=0,upper=R> num_obs_site[7];

  for (i in 1:R) {
    for (k in 1:(first[i] - 1))
      max_y[i, k] = 0;
    for (k in (last[i] + 1):7)
      max_y[i, k] = 0;
    for (k in first[i]:last[i])
      max_y[i, k] = max(y[i, 1:T, k]);
  }
  for (k in 1:7) {
    num_obs_site[k] = 0;
    for (i in 1:R)
      num_obs_site[k] = num_obs_site[k] + (y[i, 1, k] != -1);
  }
}

parameters {
  vector<upper=7>[7] alpha_lam; // Constraint for stability
  vector[7] beta;
  vector<upper=7>[R] eps_raw;   // Constraint for stability
                                // Without these constraints, some
                                // estimates become unstable maybe due
                                // to insufficient information in the
                                // data used in the BPA book.
  real<lower=0> sd_lam;
  real<lower=0> sd_p;
  vector<lower=-7,upper=7>[7] logit_p[R, T]; // Originally `lp' in the BPA book
}

transformed parameters {
  vector[R] eps;                // Abundance noise
  matrix[R, 7] log_lambda;

  eps = sd_lam * eps_raw;
  for (k in 1:7)
    for (i in 1:R)
      log_lambda[i, k] = alpha_lam[k] + eps[i];
}

model {
  // Priors
  alpha_lam ~ normal(0, sqrt(10));
  beta ~ normal(0, sqrt(10));
  eps_raw ~ normal(0, 1);

  // Weakly informative priors are used on sd_lam and sd_p,
  // instead of uniform(0, 3) used in the book
  sd_lam ~ normal(1.5, 0.75);
  sd_p ~ normal(1.5, 0.75);

  for (i in 1:R)
    for (j in 1:T)
      logit_p[i, j] ~ normal(beta, sd_p);

  // Likelihood
  for (i in 1:R)
    for (k in first[i]:last[i])
      y[i, 1:2, k] ~ bivariate_poisson_log(log_lambda[i, k],
                                           inv_logit(logit_p[i, 1:2, k]));
}

generated quantities {
  int totalN[7];            // Total pop. size across all sites
  vector[7] mean_abundance;
  vector[7] mean_N;
  vector[7] mean_detection;
  real fit = 0;
  real fit_new = 0;

  {
    int N[R, 7];          // Abundance
    real eval[R, T, 7];   // Expected values
    real y_new[R, T, 7];  // Replicate data
    vector[7] p[R, T];
    matrix[T, 7] E[R];
    matrix[T, 7] E_new[R];
    matrix[R, 7] ik_p;

    // Initialize N and ik_p
    N = rep_array(0, R, 7);
    ik_p = rep_matrix(0, R, 7);
    // Initialize E and E_new
    E[1] = rep_matrix(0, T, 7);
    E_new[1] = rep_matrix(0, T, 7);
    for (i in 2:R) {
      E[i] = E[i - 1];
      E_new[i] = E_new[i - 1];
    }
    for (i in 1:R) {
      for (j in 1:T)
        p[i, j] = inv_logit(logit_p[i, j]);

      for (k in first[i]:last[i]) {
        vector[K + 1] lp;

        for (n in 0:(max_y[i, k] - 1))
          lp[n + 1] = negative_infinity();
        for (n in max_y[i, k]:K)
          lp[n + 1] = poisson_log_lpmf(n | log_lambda[i, k])
            + binomial_lpmf(y[i, 1:T, k] | n, p[i, 1:T, k]);
        N[i, k] = categorical_rng(softmax(lp)) - 1;
        for (j in 1:T) {
          // Assess model fit using Chi-squared discrepancy
          // Compute fit statistic E for observed data
          eval[i, j, k] = p[i, j, k] * N[i, k];
          E[i, j, k] = square(y[i, j, k]
                              - eval[i, j, k]) / (eval[i, j, k] + 0.5);
          // Generate replicate data and compute fit stats for them
          y_new[i, j, k] = binomial_rng(N[i, k], p[i, j, k]);
          E_new[i, j, k] = square(y_new[i, j, k]
                                  - eval[i, j, k]) / (eval[i, j, k] + 0.5);
        }
        ik_p[i, k] = mean(p[i, 1:T, k]);
      }
    }

    for (k in 1:7) {
      totalN[k] = sum(N[, k]);
      mean_abundance[k] = mean(exp(log_lambda[, k]));
      mean_N[k] = 1.0 * totalN[k] / num_obs_site[k];
      mean_detection[k] = sum(ik_p[, k]) / num_obs_site[k];
     }
    for (i in 1:R) {
      fit = fit + sum(E[i]);
      fit_new = fit_new + sum(E_new[i]);
    }
  }
}
