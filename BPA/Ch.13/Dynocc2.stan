// Dynamic (multi-season) site-occupancy model
// This model is implemented as a hidden Markov model.

// States:
// 1: occupied
// 2: not occupied
// Observations:
// number of detection + 1
data {
  int<lower=1> nsite;                   // Number of sites
  int<lower=1> nyear;                   // Number of years
  int<lower=0,upper=2> y[nsite, nyear]; // Number of detections
}

parameters {
  real<lower=0,upper=1> psi1;
  vector<lower=0,upper=1>[nyear-1] phi;
  vector<lower=0,upper=1>[nyear-1] gamma;
  real<lower=0,upper=1> p;
}

transformed parameters {
  simplex[2] ps[2, nyear-1];    // Transition probability
  simplex[3] po[2, nyear];      // Emission probability

  for (t in 1:(nyear - 1)) {
    ps[1, t, 1] <- phi[t];
    ps[1, t, 2] <- 1.0 - phi[t];
    ps[2, t, 1] <- gamma[t];
    ps[2, t, 2] <- 1.0 - gamma[t];
  }
  for (t in 1:nyear) {
    for (r in 1:3) {
      po[1, t, r] <- exp(binomial_log(r - 1, 2, p));
      po[2, t, r] <- (r == 1);
    }
  }
}

model {
  // Priors
  // Flat priors are implicitly used on psi1, phi, gamma and p.

  // Likelihood
  // Forward algorithm derived from Stan Modeling Language
  // User's Guide  and Reference Manual
  for (i in 1:nsite) {
    real acc[2];
    vector[2] gam[nyear];

    // First year
    gam[1, 1] <- psi1 * po[1, 1, y[i, 1] + 1];
    gam[1, 2] <- (1.0 - psi1) * po[2, 1, y[i, 1] + 1];

    for (t in 2:nyear) {
      for (k in 1:2) {
        for (j in 1:2)
          acc[j] <- gam[t - 1, j] * ps[j, t - 1, k]
            * po[k, t, y[i, t] + 1];
        gam[t, k] <- sum(acc);
      }
    }
    increment_log_prob(log(sum(gam[nyear])));
  }
}

generated quantities {
  // Sample and population occupancy, growth rate and turnover
  vector[nyear] psi;                        // Occupancy probability
  int<lower=0,upper=1> z[nsite, nyear];     // Latent state of occurrence
  int<lower=0,upper=nsite> n_occ[nyear];    // Number of occupancy
  vector[nyear-1] growthr;                  // Population growth rate
  vector[nyear-1] turnover;                 // Turnover rate

  // Latent state z[,] is estimated with a full simulation
  // without using observed y[,].
  psi[1] <- psi1;
  for (k in 2:nyear)
    psi[k] <- psi[k - 1] * phi[k - 1] + (1 - psi[k - 1]) * gamma[k - 1];
  for (i in 1:nsite)
    for (k in 1:nyear)
      z[i, k] <- bernoulli_rng(psi[k]);
  for (t in 1:nyear)
    n_occ[t] <- sum(z[1:nsite, t]);
  growthr[1:(nyear - 1)] <- psi[2:nyear] ./ psi[1:(nyear - 1)];
  turnover[1:(nyear - 1)] <- (1 - psi[1:(nyear - 1)])
    .* gamma[1:(nyear - 1)] ./ psi[2:nyear];
}
