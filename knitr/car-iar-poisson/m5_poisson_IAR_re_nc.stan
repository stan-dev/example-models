data {
  int<lower=1> N_tracts; // number of census tracts in study
  int<lower=0> y[N_tracts]; // number of events per tract
  vector[N_tracts] x;  // population per tract
  int<lower=1> diag_weights[N_tracts]; // weights == num_neighbors
  int<lower=1> N_links; // number of non-zero entries in adj matrix
  int<lower=1> off_diag_coords[N_links,2]; // ij coords of off-diag entries
}
parameters {
  real beta_2;
  vector[N_tracts] re_nc;  // individual-level random effect
  real<lower=0> sigma;   // scale of random effect
  vector[N_tracts] h;  // individual-level spatial effect (IAR)
  real<lower=0> tau;  // precision param
}
transformed parameters {
  real neg_tau_div_2 = -tau * 0.5;
}
model {
  real off_diag_weight = -1.0;
  y ~ poisson_log(beta_2 * x + h + re_nc * sigma);
  beta_2 ~ normal(0, 2.5);
  tau ~ normal(0, 5);
  re_nc ~ normal(0, 1);
  sigma ~ normal(0, 5);
  for (i in 1:N_tracts) { // diagonals
    target += neg_tau_div_2 * square(h[i]) * diag_weights[i];
  }
  for (j in 1:N_links) {   // off-diagonals
    target += neg_tau_div_2
              * h[off_diag_coords[j,1]] * h[off_diag_coords[j,2]] * off_diag_weight;
  }
  target += ((N_tracts - 1) * 0.5) * log(tau);
}
generated quantities {
  vector[N_tracts] mu = exp(beta_2 * x + h);
  vector[N_tracts] SMRhat = 100.0 * mu ./ x;
}


