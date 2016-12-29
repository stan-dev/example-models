/**
 * Forward update algebra:
 * 
 * log(a' * b) = log(sum(a .* b))
 *  = log(sum(a[1] * b[1] + ... + a[N] * b[N]))
 *   = log(sum(exp(log(a[1] * b[1]) + ... + exp(log(a[N] * b[N])))))
 *   = log(sum(exp(log(a[1]) + log(b[1])) + ...)
 *   = log_sum_exp(log(a[1]) + log(b[1]), ..., log(a[N]) + log(b[N]))
 *   = log_sum_exp(log(a) + log(b))
 */
data {
  int<lower=0> N;
  vector<lower=-pi(), upper=pi()>[N] turn;
  vector<lower=0>[N] dist;
  int<lower=1> K;
}
parameters {
  simplex[K] theta[K];
  vector<lower=-pi(), upper=pi()>[K] mu;
  vector<lower=0>[K] kappa;
  positive_ordered[K] alpha;
  vector<lower=0>[K] sigma;
}  
model {
  vector[K] log_theta_tr[K];
  vector[K] lp;
  vector[K] lp_p1;

  for (k_from in 1:K)
    for (k in 1:K)
      log_theta_tr[k, k_from] = log(theta[k_from, k]);
  lp = rep_vector(-log(K), K);
  for (n in 1:N) {
    for (k in 1:K)
      lp_p1[k]
        = log_sum_exp(log_theta_tr[k] + lp)
        + von_mises_lpdf(turn[n] | mu[k], kappa[k])
        + weibull_lpdf(dist[n] | alpha[k], sigma[k]);
    lp = lp_p1;
  }
  target += log_sum_exp(lp);

  kappa ~ lognormal(0, 2);
  alpha ~ lognormal(0, 2);
  sigma ~ lognormal(0, 1);
}
