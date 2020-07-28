functions {
  vector one_comp_mm_elim_abs(real t,
                              vector y,
                              real k_a, // Dosing rate in 1/day
                              real K_m, // Michaelis-Menten constant in mg/L
                              real V_m, // Maximum elimination rate in 1/day
                              real D,
                              real V) {
    vector[1] dydt;

    real dose = 0;
    real elim = (V_m / V) * y[1] / (K_m + y[1]);

    if (t > 0)
      dose = exp(- k_a * t) * D * k_a / V;

    dydt[1] = dose - elim;

    return dydt;
  }
}

data {
  real t0;      // Initial time in days;
  vector[1] C0; // Initial concentration at t0 in mg/L

  real D;   // Total dosage in mg
  real V;   // Compartment volume in L

  int<lower=1> N_t;
  real times[N_t];   // Measurement times in days

  // Measured concentrations in effect compartment in mg/L
  real C[N_t];
}

parameters {
  real<lower=0> k_a; // Dosing rate in 1/day
  real<lower=0> K_m; // Michaelis-Menten constant in mg/L
  real<lower=0> V_m; // Maximum elimination rate in 1/day
  real<lower=0> sigma;
}

transformed parameters {
  vector[1] mu_C[N_t] = ode_bdf_tol(one_comp_mm_elim_abs, C0, t0, times,
				 1e-8, 1e-8, 1000,
				 k_a, K_m, V_m, D, V);
}

model {
  // Priors
  k_a ~ cauchy(0, 1);
  K_m ~ cauchy(0, 1);
  V_m ~ cauchy(0, 1);
  sigma ~ cauchy(0, 1);

  // Likelihood
  C ~ lognormal(log(mu_C[, 1]), sigma);
}

generated quantities {
  real C_ppc[N_t];
  for (n in 1:N_t)
    C_ppc[n] = lognormal_rng(log(mu_C[n, 1]), sigma);
}
