// Dynamic (multi-season) site-occupancy model
// This model is implemented as a hidden Markov model.

// States:
// 1: occupied
// 2: not occupied
// Observations:
// number of detection + 1
data {
  int<lower=1> nsite;                         // Number of sites
  int<lower=1> nrep;                          // Number of replicate surveys
  int<lower=1> nyear;                         // Number of years
  int<lower=0,upper=1> y[nsite, nrep, nyear]; // Detection history
}

transformed data {
  int<lower=1,upper=nrep+1> sum_y[nsite, nyear];  // sum + 1
  int ny_minus_1 = nyear - 1;

  for (i in 1:nsite)
    for (k in 1:nyear)
      sum_y[i, k] = sum(y[i, 1:nrep, k]) + 1;
}

parameters {
  real<lower=0,upper=1> psi1;
  vector<lower=0,upper=1>[nyear] p;
  simplex[2] ps[2, nyear-1];     // Transition probability
  // This is equivalent to the following.
  //  ps[1, t, 1] = phi[t];
  //  ps[1, t, 2] = 1.0 - phi[t];
  //  ps[2, t, 1] = gamma[t];
  //  ps[2, t, 2] = 1.0 - gamma[t];
}

transformed parameters {
  simplex[nrep + 1] po[2, nyear];  // Log emission probability

  for (t in 1:nyear) {
    for (r in 1:(nrep + 1)) {
      po[1, t, r] = exp(binomial_lpmf(r - 1 | nrep, p[t]));
      po[2, t, r] = (r == 1);
    }
  }
}

model {
  // Priors
  // Flat priros Uniform(0, 1) are implicitly used on psi1, p and ps.

  // Likelihood
  // This implementation of the forward algorithm is derived from
  // Stan Modeling Language User's Guide and Reference Manual.
  for (i in 1:nsite) {
    real acc[2];
    vector[2] gam[nyear];

    // First year
    gam[1, 1] = psi1 * po[1, 1, sum_y[i, 1]];
    gam[1, 2] = (1 - psi1) * po[2, 1, sum_y[i, 1]];

    for (t in 2:nyear) {
      for (k in 1:2) {
        for (j in 1:2)
          acc[j] = gam[t - 1, j] * ps[j, t - 1, k] * po[k, t, sum_y[i, t]];
        gam[t, k] = sum(acc);
      }
    }
    target += log(sum(gam[nyear]));
  }
}

generated quantities {
  // Sample and population occupancy, growth rate and turnover
  vector<lower=0,upper=1>[nyear] psi;       // Occupancy probability
  vector<lower=0,upper=1>[ny_minus_1] phi;   // Survival probability
  vector<lower=0,upper=1>[ny_minus_1] gamma; // Colonization probability
  int<lower=0,upper=1> z[nsite, nyear];     // Latent state of occurrence
  int<lower=0,upper=nsite> n_occ[nyear];    // Number of occupancy
  vector[nyear-1] growthr;                  // Population growth rate
  vector[nyear-1] turnover;                 // Turnover rate

  for (k in 1:ny_minus_1) {
    phi[k] = ps[1, k, 1];
    gamma[k] = ps[2, k, 1];
  }
  psi[1] = psi1;
  for (k in 2:nyear)
    psi[k] = psi[k - 1] * phi[k - 1] + (1 - psi[k - 1]) * gamma[k - 1];
  for (i in 1:nsite)
    for (k in 1:nyear)
      z[i, k] = bernoulli_rng(psi[k]);
  for (t in 1:nyear)
    n_occ[t] = sum(z[1:nsite, t]);
  growthr[:ny_minus_1] = psi[2:] ./ psi[:ny_minus_1];
  turnover[:ny_minus_1] = (1 - psi[:ny_minus_1]) .* gamma[:ny_minus_1] ./ psi[2:];
}
