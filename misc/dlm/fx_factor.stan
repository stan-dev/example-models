// Kalman Filter (Multivariate Form)
// - no missing values
data {
  // system matrices
  int r;
  int T;
  matrix[r, T] y;
  vector[1] m0;
  cov_matrix[1] C0;
}
parameters {
  matrix[1, 1] G;
  vector[r - 1] lambda;
  vector<lower=0.0>[r] V;
  cov_matrix[1] W;
}
transformed parameters {
  matrix[1, r] F;
  F[1, 1] = 1;
  for (i in 1:(r - 1)) {
    F[1, i + 1] = lambda[i];
  }
}
model {
  matrix[1,1] identity;
  identity <- diag_matrix(rep_vector(1.0,1)); 
  W ~ inv_wishart(1, identity);
  to_vector(G) ~ normal(0, 10);
  lambda ~ normal(0, 10);
  V ~ normal(0, 10);
  y ~ gaussian_dlm_obs(F, G, V, W, m0, C0);
}