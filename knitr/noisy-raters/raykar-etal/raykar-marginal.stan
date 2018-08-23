/**
 * Implementation of joint estimation of logistic regression and
 * annotator sensitivity and specifity.
 *
 * TODO:  hierarchical priors for alpha, beta
 * TODO:  precompute counts for sum(y[ , r] == 1) and multiply
 *        in the N/R loops and the eltwise multiply and sum
 */
data {
  int<lower=0> N;                 // # items (reviews)
  int<lower=0> R;                 // # raters (heuristics)
  int<lower=0> D;                 // # of predictors for each item
  matrix[N, D] x;                 // [n, d] predictors d for item n
  int<lower=0, upper=1> y[N, R];  // outcome: 1 success, 0 failure
}
parameters {
  vector[D] w;                         // logistic regression coeffs
  real w0;                             // intercept
  vector<lower=0, upper=1>[R] alpha;   // sensitivity
  vector<lower=0, upper=1>[R] beta;    // specificity
}
model {
  vector[N] logit_z_hat = w0 + x * w;
  vector[N] log_z_hat = log_inv_logit(logit_z_hat);
  vector[N] log1m_z_hat = log1m_inv_logit(logit_z_hat);

  vector[R] log_alpha = log(alpha);
  vector[R] log1m_alpha = log1m(alpha);
  vector[R] log_beta = log(beta);
  vector[R] log1m_beta = log1m(beta);

  // prior
  w ~ normal(0, 2);
  w0 ~ normal(0, 5);
  alpha ~ beta(10, 1);
  beta ~ beta(10, 1);

  // likelihood
  for (n in 1:N) {
    real pos_sum = log_z_hat[n];
    real neg_sum = log1m_z_hat[n];
    for (r in 1:R) {
      if (y[n, r] == 1) {
        pos_sum += log_alpha[r];
        neg_sum += log1m_beta[r];
      } else {
        pos_sum += log1m_alpha[r];
        neg_sum += log_beta[r];
      }
    }
    target += log_sum_exp(pos_sum, neg_sum);
  }
}
generated quantities {
  vector[N] Pr_z_eq_1;
  {
    vector[N] logit_z_hat = w0 + x * w;
    vector[N] log_z_hat = log_inv_logit(logit_z_hat);
    vector[N] log1m_z_hat = log1m_inv_logit(logit_z_hat);

    vector[R] log_alpha = log(alpha);
    vector[R] log1m_alpha = log1m(alpha);
    vector[R] log_beta = log(beta);
    vector[R] log1m_beta = log1m(beta);

    for (n in 1:N) {
      real pos_sum = log_z_hat[n];
      real neg_sum = log1m_z_hat[n];
      for (r in 1:R) {
        if (y[n, r] == 1) {
          pos_sum += log_alpha[r];
          neg_sum += log1m_beta[r];
        } else {
          pos_sum += log1m_alpha[r];
          neg_sum += log_beta[r];
        }
      }
      Pr_z_eq_1[n] = softmax([pos_sum, neg_sum]')[1];
    }
  }
}
