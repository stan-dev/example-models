data {
  int N; // Number of observations
  array[N] real y;
}
parameters {
  // Parameters of measurement model
  ordered[3] mu;
  real<lower=0.0> sigma;
  
  // Initial state
  simplex[3] rho;
  
  // Rows of the transition matrix
  simplex[2] t1;
  simplex[3] t2;
  simplex[2] t3;
}
transformed parameters {
  matrix[3, 3] Gamma = rep_matrix(0, 3, 3);
  matrix[3, N] log_omega;
  
  // Build the transition matrix
  Gamma[1, 1 : 2] = t1';
  Gamma[2,  : ] = t2';
  Gamma[3, 2 : 3] = t3';
  
  // Compute the log likelihoods in each possible state
  for (n in 1 : N) {
    // The observation model could change with n, or vary in a number of
    //  different ways (which is why log_omega is passed in as an argument)
    log_omega[1, n] = normal_lpdf(y[n] | mu[1], sigma);
    log_omega[2, n] = normal_lpdf(y[n] | mu[2], sigma);
    log_omega[3, n] = normal_lpdf(y[n] | mu[3], sigma);
  }
}
model {
  mu ~ normal(0, 10);
  sigma ~ normal(0, 1);
  
  rho ~ dirichlet([10, 1, 1]);
  
  t1 ~ dirichlet([1, 1]);
  t2 ~ dirichlet([1, 1, 1]);
  t3 ~ dirichlet([1, 1]);
  
  target += hmm_marginal(log_omega, Gamma, rho);
}
generated quantities {
  matrix[3, N] hidden_probs = hmm_hidden_state_prob(log_omega, Gamma, rho);
  array[N] int y_sim = hmm_latent_rng(log_omega, Gamma, rho);
}
