data {
  int n;
  matrix[1, n] y;
  vector[1] m0;
  matrix[1, 1] C0;
}
transformed data {
  matrix[1, 1] F;
  matrix[1, 1] G;
  F = rep_matrix(1, 1, 1);
  G = rep_matrix(1, 1, 1);
}
parameters {
  real<lower = 0> sigma_y;
  real<lower = 0> sigma_theta;
}
model {
  matrix[1, 1] V;
  matrix[1, 1] W;
  V[1, 1] =  pow(sigma_y, 2);
  W[1, 1] = pow(sigma_theta, 2);
  y ~ gaussian_dlm_obs(F, G, V, W, m0, C0);
}
