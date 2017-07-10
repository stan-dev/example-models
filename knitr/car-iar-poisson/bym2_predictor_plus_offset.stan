data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> y[N];              // count outcomes
  vector[N] x;                    // predictor
  vector<lower=0>[N] E;           // exposure
}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  real beta0;                // intercept
  real beta1;                // slope

  real<lower=0> sigma_random;   // precision of heterogeneous effects
  real<lower=0,upper=1> mixing_parameter;     // precision of spatial effects

  vector[N] theta_std;       // standardized heterogeneous effects
  vector[N - 1] phi_std_raw; // raw, standardized spatial effects
}
transformed parameters {

  vector[N] phi;
  vector[N] random;

  phi[1:(N - 1)] = phi_std_raw;
  phi[N] = -sum(phi_std_raw);

  // non-centered parameterisation
  random =  sigma_random *(sqrt(mixing_parameter)*theta_std + sqrt(1-mixing_parameter)*phi);
}
model {
  y ~ poisson_log(log_E + beta0 + beta1 * x + random);

  target += -0.5 * dot_self(phi[node1] - phi[node2]);

  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  theta_std ~ normal(0, 1);
  sigma_random ~ normal(0,5);
  mixing_parameter ~ beta(0.5,0.5);
}
generated quantities {
  vector[N] mu = exp(log_E + beta0 + beta1 * x +random);
  //real psi = sd(phi) / (sd(theta) + sd(phi));  // proportion spatial variation
}
