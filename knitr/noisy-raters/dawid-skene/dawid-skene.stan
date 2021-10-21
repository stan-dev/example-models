/*
 * Implementation of Dawid and Skene's model for measurement error
 * among categoical raters.  The model is extended for Bayes fixed priors.
 * Discrete parameters for category of items marginalized out with their
 * expectation returned on the log scale in derived quantity log_q_z.
 *
 * Dawid, A. P., & Skene, A. M. (1979). Maximum likelihood estimation
 * of observer error-rates using the EM algorithm. Applied statistics,
 * 20--28.
 */
data {
  int<lower=2> K; // number of categories
  int<lower=1> I; // number of items
  int<lower=1> J; // number of coders
  array[I, J] int<lower=1, upper=K> y; // label for observation n
  vector<lower=0>[K] alpha; // prior for prevalence (positive)
  array[K] vector<lower=0>[K] beta; // prior for coder responses (positive)
}
parameters {
  simplex[K] pi; // prevalence of categories
  array[J, K] simplex[K] theta; // response of anotator j to category k
}
model {
  // log prior:  log p(theta, pi)
  pi ~ dirichlet(alpha);
  for (j in 1 : J) {
    for (k in 1 : K) {
      theta[j, k] ~ dirichlet(beta[k]);
    }
  }
  
  // log marginal likelihood: log p(y | theta, pi)
  for (i in 1 : I) {
    vector[K] log_q = log(pi);
    for (j in 1 : J) {
      log_q += to_vector(log(theta[j,  : , y[i, j]]));
    }
    target += log_sum_exp(log_q);
  }
}
generated quantities {
  // normalized log posterior: log Pr[z | y]
  array[I] vector[K] log_Pr_z;
  for (i in 1 : I) {
    vector[K] log_q = log(pi);
    for (j in 1 : J) {
      log_q += to_vector(log(theta[j,  : , y[i, j]]));
    }
    log_Pr_z[i] = log_q - log_sum_exp(log_q);
  }
}
