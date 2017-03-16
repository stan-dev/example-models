data {
  int<lower=1> K;
  int<lower=1> D;
  int<lower=0> N;
  int<lower=0,upper=1> y[N,D];
  vector[K] x[N];
}
parameters {
  matrix[D,K] beta;
  cholesky_factor_corr[D] L_Omega;
  real<lower=0,upper=1> u[N,D]; // nuisance that absorbs inequality constraints
}
model {
  L_Omega ~ lkj_corr_cholesky(4);
  to_vector(beta) ~ normal(0, 5);
  // implicit: u is iid standard uniform a priori
  { // likelihood
    for (n in 1:N) {
      vector[D] mu;
      vector[D] z;
      real prev;
      mu = beta * x[n];
      prev = 0;
      for (d in 1:D) { // Phi and inv_Phi may overflow and / or be numerically inaccurate
        real bound; // threshold at which utility = 0
        bound = Phi( -(mu[d] + prev) / L_Omega[d,d]  );
        if (y[n,d] == 1) {
          real t;
          t = bound + (1 - bound) * u[n,d];
          z[d] = inv_Phi(t);       // implies utility is positive
          target += log1m(bound);  // Jacobian adjustment
        }
        else {
          real t;
          t = bound * u[n,d];
          z[d] = inv_Phi(t);     // implies utility is negative
          target += log(bound);  // Jacobian adjustment
        }
        if (d < D) prev = L_Omega[d+1,1:d] * head(z, d);
        // Jacobian adjustments imply z is truncated standard normal
        // thus utility --- mu + L_Omega * z --- is truncated multivariate normal
      }
    }
  }
}
generated quantities {
  corr_matrix[D] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);  
}
