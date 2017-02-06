// Open-population binomial-mixuture model

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
   * @param p          Detection probability
   *
   * return Log probability
  */
  real bivariate_poisson_log_lpmf(int[] n, real log_lambda, real p) {
    real s[min(n) + 1];
    real log_theta_1 = log_lambda + log(p) + log1m(p);
    real log_theta_0 = log_lambda + log(p) * 2;

    if (size(n) != 2)
      reject("Size of n must be 2.");
    if (p < 0 || p > 1)
      reject("p must be in [0,1].");
    for (u in 0:min(n))
      s[u + 1] = poisson_log_lpmf(n[1] - u | log_theta_1)
               + poisson_log_lpmf(n[2] - u | log_theta_1)
               + poisson_log_lpmf(u | log_theta_0);
    return log_sum_exp(s);
  }
}

data {
  int<lower=1> R;                // Number of sites
  int<lower=1> T;                // Number of replications; fixed as 2
  int<lower=-1> y[R, 2, 7];      // Counts (-1:NA)
  int<lower=1,upper=7> first[R]; // First occasion
  int<lower=1,upper=7> last[R];  // Last occasion
  int<lower=1> K;
}

transformed data {
  int<lower=0> max_y[R, 7];

  for (i in 1:R) {
    for (k in 1:(first[i] - 1))
      max_y[i, k] = 0;
    for (k in (last[i] + 1):7)
      max_y[i, k] = 0;
    for (k in first[i]:last[i])
      max_y[i, k] = max(y[i, 1:T, k]);
  }
}

parameters {
  vector[7] alpha_lam;
  vector<lower=0,upper=1>[7] p;  // Capture probability
}

model {
  // Priors
  alpha_lam ~ normal(0, 10);
  // A flat prior Uniform(0, 1) is implicitly used on p.

  // Likelihood
  for (i in 1:R)
    for (k in first[i]:last[i])
      y[i, 1:T, k] ~ bivariate_poisson_log(alpha_lam[k], p[k]);
}

generated quantities {
  int totalN[7];
  real fit = 0;
  real fit_new = 0;
  vector[7] mean_abundance;

  {
    int N[R, 7];            // Abundance
    real eval[R, 7];        // Expected values
    int y_new[R, T, 7];
    matrix[T, 7] E[R];
    matrix[T, 7] E_new[R];

    // Initialize N, E and E_new
    N = rep_array(0, R, 7);
    E[1] = rep_matrix(0, T, 7);
    E_new[1] = rep_matrix(0, T, 7);
    for (i in 2:R) {
      E[i] = E[i - 1];
      E_new[i] = E_new[i - 1];
    }
    for (i in 1:R) {
      for (k in first[i]:last[i]) {
        vector[K + 1] lp;

        for (n in 0:(max_y[i, k] - 1))
          lp[n + 1] = negative_infinity();
        for (n in max_y[i, k]:K)
          lp[n + 1] = poisson_log_lpmf(n | alpha_lam[k])
            + binomial_lpmf(y[i, 1:T, k] | n, p[k]);
        N[i, k] = categorical_rng(softmax(lp)) - 1;
        eval[i, k] = p[k] * N[i, k];
        for (j in 1:T) {
          // Assess model fit using Chi-squared discrepancy
          // Compute fit statistic E for observed data
          E[i, j, k] = square(y[i, j, k] - eval[i, k]) / (eval[i, k] + 0.5);
          // Generate replicate data and
          // Compute fit statistic E_new for replicate data
          y_new[i, j, k] = binomial_rng(N[i, k], p[k]);
          E_new[i, j, k] = square(y_new[i, j, k] - eval[i, k]) / (eval[i, k] + 0.5);
        }
      }
    }
    for (k in 1:7)
      totalN[k] = sum(N[1:R, k]);  // Total pop. size across all sites
    for (i in 1:R) {
      fit = fit + sum(E[i]);
      fit_new = fit_new + sum(E_new[i]);
    }
  }
  mean_abundance = exp(alpha_lam);
}
