data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
  vector<lower=0>[N] E;           // exposure
}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  real<lower=0> tau_phi;     // precision of spatial effects
  vector[N - 1] phi_std_raw; // raw, standardized spatial effects
}
transformed parameters {
  real<lower=0> sigma_phi = inv(sqrt(tau_phi));      // convert precision to sigma
  vector[N] phi;
  phi[1:(N - 1)] = phi_std_raw;
  phi[N] = -sum(phi_std_raw);
  phi = phi * sigma_phi;    // non-centered parameterization
}
model {
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  tau_phi ~ gamma(1, 1);            // Carlin WinBUGS priors
}
generated quantities {
  real beta0;                // intercept
  real beta1;                // slope
  real<lower=0> tau_theta;   // precision of heterogeneous effects
  real<lower=0> sigma_theta;
  vector[N] theta;
  vector[N] theta_std;

  vector[N] x;
  int y[N];
  
  beta0 = normal_rng(0, 1);
  beta1 = normal_rng(0, 1);
  for (i in 1:N) {
    theta_std[i] = normal_rng(0, 1);
  }
  tau_theta = gamma_rng(3.2761, 1.81);  // Carlin WinBUGS priors
  sigma_theta = inv(sqrt(tau_theta));  // convert precision to sigma
  theta = theta_std * sigma_theta;         // non-centered parameterization

  for (i in 1:N) {
    x[i] = normal_rng(0, 1);
    y[i] = poisson_log_rng(log_E[i] + beta0 + beta1 * x[i] + phi[i] + theta[i]);
  }
}
