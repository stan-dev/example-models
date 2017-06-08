data {
  int<lower=1> N_areas;
  int<lower=0> y[N_areas]; // number of events per area
  vector[N_areas] x_aff;  // covariate: pct pop employed in AFF
  vector[N_areas] e_pop;  // exposure: population
  int<lower=1> diag_weights[N_areas]; // weights == num_neighbors
  int<lower=1> N_links; // number of non-zero entries in adj matrix
  int<lower=1> off_diag_coords[N_links,2]; // ij coords of off-diag entries
}
transformed data {
  vector[N_areas] log_e_pop = log(e_pop);
}
parameters {
  real beta_2;   // slope
  vector[N_areas] h;  // individual-level spatial effect (IAR)
  real<lower=0> tau;  // precision param
  vector[N_areas] re_nc;  // individual-level random effect
  real<lower=0> sigma;   // scale of random effect
}
transformed parameters {
  real neg_tau_div_2 = -tau * 0.5;
}
model {
  real off_diag_weight = -1.0;
  y ~ poisson_log(beta_2 * x_aff + log_e_pop + h + re_nc * sigma);
  beta_2 ~ normal(0, 1);
  tau ~ gamma(1.0, 1.0);  // following Banerjee et al 2003
  re_nc ~ normal(0, 1);
  sigma ~ gamma(3.3, 1.8);   // following Banerjee et al 2003
  target += ((N_areas - 1) * 0.5) * log(tau);
  for (i in 1:N_areas) { // diagonals
    target += neg_tau_div_2 * square(h[i]) * diag_weights[i];
  }
  for (j in 1:N_links) {   // off-diagonals
    target += neg_tau_div_2
              * h[off_diag_coords[j,1]] * h[off_diag_coords[j,2]] * off_diag_weight;
  }
}
generated quantities {
  real beta_1;
  vector[N_areas] eta = h + re_nc * sigma;
  vector[N_areas] mu = exp(beta_2 * x_aff + log_e_pop + eta);
  vector[N_areas] intercepts;
  for (i in 1:N_areas) {
    intercepts[i] = y[i] - mu[i];
  }
  beta_1 = mean(intercepts);
}
