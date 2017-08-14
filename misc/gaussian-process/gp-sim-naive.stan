// Sample from Gaussian process with a manually-built covariance matrix K
// parameterized with an exponentiated quadratic kernel.  
// Fixed kernel hyperparameters: rho=1, alpha=1, sigma=sqrt(0.1)

data {
  int<lower=1> N;
  real x[N];
}
transformed data {
  matrix[N, N] K;
  vector[N] mu = rep_vector(0, N);
  for (i in 1:(N - 1)) {
    K[i, i] = 1 + 0.1;
    for (j in (i + 1):N) {
      K[i, j] = exp(-0.5 * pow(x[i] - x[j],2));
      K[j, i] = K[i, j];
    }
  }
  K[N, N] = 1 + 0.1;
}
parameters {
  vector[N] y;
}
model {
  y ~ multi_normal(mu, K);
}
