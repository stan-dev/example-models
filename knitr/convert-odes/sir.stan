// Simple SIR model inspired by the presentation in
// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3380087/pdf/nihms372789.pdf

functions {
  vector simple_SIR(real t,
                    vector y,
                    real beta,    // water contact rate
                    real kappa,   // C_{50}
                    real gamma,   // recovery rate
                    real xi,      // bacteria production rate
                    real delta) { // bacteria removal rate
    vector[4] dydt;

    dydt[1] = -beta * y[4] / (y[4] + kappa) * y[1];
    dydt[2] = beta * y[4] / (y[4] + kappa) * y[1] - gamma * y[2];
    dydt[3] = gamma * y[2];
    dydt[4] = xi * y[2] - delta * y[4];

    return dydt;
  }
}

data {
  int<lower=0> N_t;
  real t[N_t];
  vector[4] y0;
  int stoi[N_t];
  real B[N_t];
}

transformed data {
  real t0 = 0;
  real<lower=0> kappa = 1000000;
}

parameters {
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> xi;
  real<lower=0> delta;
}

transformed parameters {
  vector<lower=0>[4] y[N_t] = ode_rk45_tol(simple_SIR, y0, t0, t,
					   1e-6, 1e-6, 1000,
					   beta, kappa, gamma, xi, delta);
}
  
model {
  vector[N_t] y_diff;
  y_diff[1] = y0[1] - y[1, 1];
  for (n in 2:N_t)
    y_diff[n] = y[n - 1, 1] - y[n, 1];
  
  beta ~ cauchy(0, 2.5);
  gamma ~ cauchy(0, 1);
  xi ~ cauchy(0, 25);
  delta ~ cauchy(0, 1);

  stoi ~ poisson(y_diff);
  B ~ lognormal(log(y[, 4]), 0.15);  
}
