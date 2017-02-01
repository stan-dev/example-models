// Zero-inflated Poisson binomial-mixture model

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
   * @param n      Number of observed individuals
   * @param lambda Poisson mean of population size
   * @param p      Detection probability
   *
   * return Log probability
   */
  real bivariate_poisson_lpmf(int[] n, real lambda, real p) {
    real s[min(n) + 1];
    real theta_1 = lambda * p * (1 - p);
    real theta_0 = lambda * p * p;

    if (lambda < 0) {
      reject("lambda must be non-negative.");
    } else if (p < 0 || p > 1) {
      reject("p must be in [0,1].");
    }
    for (u in 0:min(n))
      s[u + 1] = poisson_lpmf(n[1] - u | theta_1)
               + poisson_lpmf(n[2] - u | theta_1)
               + poisson_lpmf(u | theta_0);
    return log_sum_exp(s);
  }

  /**
   * Return log probability of Zero-inflated Poisson Binomial Mixture
   *
   * @param y          Count
   * @param max_y_site Max. number of counts on the site
   * @param N          Population size
   * @param oemga      Inclusion probability
   * @param p          Detection probability
   *
   * @return Log probability
   */
  real zipbin_lpmf(int[] y, int max_y_site,
                   int n, real omega, real log_lambda, real p) {
    real lp;

    if (max_y_site) {
      if (max(y) > n) {
        lp = negative_infinity();
      } else {
        lp = bernoulli_lpmf(1 | omega)
          + poisson_log_lpmf(n | log_lambda)
          + binomial_lpmf(y | n, p);
      }
    } else {
      lp = log_sum_exp(bernoulli_lpmf(0 | omega),
                       bernoulli_lpmf(1 | omega)
                       + poisson_log_lpmf(n | log_lambda)
                       + binomial_lpmf(y | n, p));
    }
    return lp;
  }
}

data {
  int<lower=1> R;                // Number of sites
  // int<lower=1> T;                // Number of replications; fixed as 2
  int<lower=-1> y[R, 2, 7];      // Counts (-1:NA)
  int<lower=1,upper=7> first[R]; // First occasion
  int<lower=1,upper=7> last[R];  // Last occasion
  int<lower=0> K;                // Upper bounds of population size
}

transformed data {
  int<lower=0> max_y[R, 7];
  int<lower=0> max_y_site[R];
  int T = 2;

  for (i in 1:R) {
    for (k in 1:(first[i] - 1))
      max_y[i, k] = 0;
    for (k in (last[i] +1 ):7)
      max_y[i, k] = 0;
    for (k in first[i]:last[i])
      max_y[i, k] = max(y[i, 1:T, k]);
    max_y_site[i] = max(max_y[i]);
  }
}

parameters {
  real<lower=0,upper=1> omega;   // Suitability
  vector[7] alpha_lam;           // Log abundance
  vector<lower=0,upper=1>[7] p;  // Captue probability
}

model {
  // Priors
  // Implicit flat priors [0, 1] are used on omega and p.
  alpha_lam ~ normal(0, 10);

  // Likelihood
  for (i in 1:R) {
    if (max_y_site[i]) {
      real lp = bernoulli_lpmf(1 | omega);

      for (k in first[i]:last[i])
        lp = lp + bivariate_poisson_lpmf(y[i, 1:2, k] |
                                         exp(alpha_lam[k]), p[k]);
      target += lp;
    } else {
      real lp[2];

      lp[1] = bernoulli_lpmf(0 | omega);
      lp[2] = bernoulli_lpmf(1 | omega);
      for (k in first[i]:last[i])
        lp[2] = lp[2] + bivariate_poisson_lpmf(y[i, 1:2, k] |
                                               exp(alpha_lam[k]), p[k]);
      target += log_sum_exp(lp);
    }
  }
}

generated quantities {
  int totalN[7];         // Total pop. size across all sites
  real fit = 0;
  real fit_new = 0;
  vector[7] mean_abundance;

  {
    int N[R, 7];         // Latent abundance state
    real eval[R, 7];     // Expected values
    int y_new[R, T, 7];
    matrix[T, 7] E[R];
    matrix[T, 7] E_new[R];

    N = rep_array(0, R, 7);
    for (i in 1:R) {
      real p_unobs;  // Prob. site is suitable but no indiv. observed.

      E[i] = rep_matrix(0, T, 7);
      E_new[i] = rep_matrix(0, T, 7);

      for (k in first[i]:last[i]) {
        vector[K + 1] lp;

        for (n in 0:K)
          lp[n + 1] = zipbin_lpmf(y[i, 1:T, k] | max_y_site[i],
                                    n, omega, alpha_lam[k], p[k]);
        N[i, k] = categorical_rng(softmax(lp)) - 1;
      }

      if (max_y_site[i] == 0) {  // Unobserved
        p_unobs = omega * exp(binomial_lpmf(0 | N[i], p))^T;
        if (bernoulli_rng(p_unobs) == 0) {
          // Site is not suitable
          for (k in first[i]:last[i])
            N[i, k] = 0;
        }
      }

      for (k in first[i]:last[i]) {
        eval[i, k] = p[k] * N[i, k];
        for (j in 1:T) {
          // Assess model fit using Chi-squared discrepancy
          // Compute fit statistic E for observed data
          E[i, j, k] = square(y[i, j, k] - eval[i, k])
            / (eval[i, k] + 0.5);
          // Generate replicate data and compute fit stats for them
          y_new[i, j, k] = binomial_rng(N[i, k], p[k]);
          E_new[i, j, k] = square(y_new[i, j, k] - eval[i, k])
            / (eval[i, k] + 0.5);
        }
      }
    }
    for (k in 1:7)
      totalN[k] = sum(N[1:R, k]);
    for (i in 1:R) {
      fit = fit + sum(E[i]);
      fit_new = fit_new + sum(E_new[i]);
    }
  }
  mean_abundance = exp(alpha_lam);
}
