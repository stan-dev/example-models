data {
  int<lower=1> N_tracts; // number of census tracts in study
  int<lower=0> y[N_tracts];
  vector[N_tracts] x;
  int<lower=1> D_W[N_tracts]; // weights == num_neighbors
  int<lower=1> N_links; // number of non-zero entries in adj matrix
  int<lower=1> off_diag_coords[N_links,2]; // ij coords of off-diag entries
}
transformed data {
  vector[N_tracts] x_nonzero = x + 0.001;
}
parameters {
  real beta_2;
  vector[N_tracts] h;  // CAR component
  real<lower=0> tau;  // precision param
}
transformed parameters {
  real neg_tau_div_2 = -tau / 2;
}
model {
  y ~ poisson_log(h + beta_2 * x_nonzero);
  beta_2 ~ normal(0, 2.5);
  tau ~ normal(0, 5);
  for (k in 1:N_links) {   // off-diagonals
    target += neg_tau_div_2 *
      (-1.0) * h[off_diag_coords[k,1]] * h[off_diag_coords[k,2]];
  }
  for (k in 1:N_tracts) { // diagonals
    target += neg_tau_div_2 * square(h[k]) * D_W[k];
  }
  target += ((N_tracts - 1) / 2.0) * log(tau);
}
generated quantities {
  vector[N_tracts] mu = exp(h + beta_2 * x_nonzero);
}


