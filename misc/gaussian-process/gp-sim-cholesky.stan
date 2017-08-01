// Sample from Gaussian process
// Fixed covar function: eta_sq=1, rho_sq=1, sigma_sq=0.1

data {
  int<lower=1> N;
  int<lower=1> D;
  vector[D] x[N];
}
transformed data {
  vector[N] mu = rep_vector(0, N);
  matrix[N, N] L;
  {
    matrix[N, N] K = cov_exp_quad(x, 1, 1);
    for (i in 1:N) 
      K[i, i] = K[i, i] + 0.1;
    L = cholesky_decompose(K);
  }
}
parameters {
  vector[N] z;
}
model {
  z ~ normal(0, 1);
}
generated quantities {
  vector[N] y;
  y = mu + L * z;
}
