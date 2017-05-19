data {
  int<lower=1> N_tracts; // number of census tracts in study
  int<lower=0> y[N_tracts];
  vector[N_tracts] x;
  int<lower=1> diag_weights[N_tracts]; // weights == num_neighbors
  int<lower=1> N_links; // number of non-zero entries in adj matrix
  int<lower=1> off_diag_coords[N_links,2]; // ij coords of off-diag entries
}
parameters {
  real beta_2;
  vector[N_tracts] h;  // individual-level spatial effect (IAR)
  real<lower=0> tau;  // precision param
}
transformed parameters {
  real neg_tau_div_2 = -tau / 2;
}
model {
  real off_diag_weight = -1.0;
  y ~ poisson_log(beta_2 * x + h);
  beta_2 ~ normal(0, 2.5);
  tau ~ normal(0, 5);
  for (i in 1:N_tracts) { // diagonals
    target += neg_tau_div_2 * square(h[i]) * diag_weights[i];
  }
  for (j in 1:N_links) {   // off-diagonals
    target += neg_tau_div_2 *
    h[off_diag_coords[j,1]] * h[off_diag_coords[j,2]] * off_diag_weight;
  }
  target += ((N_tracts - 1) / 2.0) * log(tau);
}
generated quantities {
  vector[N_tracts] mu = exp(beta_2 * x + h);
}


