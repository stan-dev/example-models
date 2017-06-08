data {
  int<lower=1> N_areas;
  int<lower=0> y[N_areas]; // number of events per area
  vector[N_areas] x;  // population per area
  int<lower=1> diag_weights[N_areas]; // weights == num_neighbors
  int<lower=1> N_links; // number of non-zero entries in adj matrix
  int<lower=1> off_diag_coords[N_links,2]; // ij coords of off-diag entries
}
parameters {
  real beta_2;
  vector[N_areas] re_nc;  // individual-level random effect
  real<lower=0> sigma;   // scale of random effect
  vector[N_areas] h;  // individual-level spatial effect (IAR)
  real<lower=0> tau;  // precision param
}
transformed parameters {
  real neg_tau_div_2 = -tau * 0.5;
}
model {
  real off_diag_weight = -1.0;
  y ~ poisson_log(beta_2 * x + h + re_nc * sigma);
  beta_2 ~ normal(0, 2.5);
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
  vector[N_areas] beta_1s;
  vector[N_areas] mu = exp(beta_2 * x + h + re_nc * sigma);
  for (i in 1:N_areas) {
    beta_1s[i] = y[i] - mu[i];
  }
  beta_1 = mean(beta_1s);
}


