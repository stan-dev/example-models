// K = 2 fixed
data {
  int<lower=0> N;
  vector<lower=-pi(), upper=pi()>[N] turn;
  vector<lower=0>[N] dist;
}
parameters {
  vector<lower = 0, upper = 1> theta[2];
  vector<lower=-pi(), upper=pi()>[2] mu;
  vector<lower=0>[2] kappa;
  positive_ordered[2] alpha;
  vector<loweer=0>[2] sigma;
}  
model {
  vector[2] log_theta;
  vector[2] log1m_theta;
  vector[2] lp;

  log_theta = log(theta);
  for (k in 1:2) log1m_theta[2] = log1m(theta);
  lp = rep_vector(-log(2), 2);
  for (n in 1:N) {
    for (k in 1:2)
      lp[k]
        = log_sum_exp((k == 1 ? log_theta[1] : log1m_theta[1])+ lp)
        + von_mises_lpdf(turn[n] | mu[k], kappa[k])
        + lognormal_lpdf(dist[n] | alpha[k], sigma[k]);
  }
  target += log_sum_exp(lp);

  kappa ~ lognormal(0, 2);
  alpha ~ lognormal(0, 2);
  sigma ~ lognormal(0, 1);
}
