// Jolly-Seber model as s multistate model

//--------------------------------------
// States (S):
// 1 not yet entered
// 2 alive
// 3 dead
// Observations (O):
// 1 seen
// 2 not seen
//--------------------------------------

data {
  int<lower=0> M;                               // Augmented sample size
  int<lower=0> n_occasions;                     // Number of capture occasions
  int<lower=1,upper=2> y[M, n_occasions];       // Augmented capture-history
}

parameters {
  vector<lower=0,upper=1>[n_occasions-1] gamma; // Removal entry probabilities
  real<lower=0,upper=1> mean_phi;               // Mean survival
  real<lower=0,upper=1> mean_p;                 // Mean capture
}

transformed parameters {
  vector<lower=0,upper=1>[n_occasions-1] phi; // Survival probability
  vector<lower=0,upper=1>[n_occasions-1] p;   // Capture probability
  simplex[3] ps[3, M, n_occasions-1];
  simplex[2] po[3, M, n_occasions-1];

  // Priors and constraints
  for (t in 1:(n_occasions - 1)) {
    phi[t] <- mean_phi;
    p[t] <- mean_p;
  }

  // Define state-transition and observation matrices
  for (i in 1:M) {
    // Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n_occasions - 1)) {
      ps[1, i, t, 1] <- 1.0 - gamma[t];
      ps[1, i, t, 2] <- gamma[t];
      ps[1, i, t, 3] <- 0.0;
      ps[2, i, t, 1] <- 0.0;
      ps[2, i, t, 2] <- phi[t];
      ps[2, i, t, 3] <- 1 - phi[t];
      ps[3, i, t, 1] <- 0.0;
      ps[3, i, t, 2] <- 0.0;
      ps[3, i, t, 3] <- 1.0;

      // Define probabilities of O(t) given S(t)
      po[1, i, t, 1] <- 0.0;
      po[1, i, t, 2] <- 1.0;
      po[2, i, t, 1] <- p[t];
      po[2, i, t, 2] <- 1.0 - p[t];
      po[3, i, t, 1] <- 0.0;
      po[3, i, t, 2] <- 1.0;
    } //t
  } //i
}

model {
  real acc[3];
  vector[3] gam[n_occasions];

  // Priors
  gamma ~ uniform(0, 1);
  mean_phi ~ uniform(0, 1);
  mean_p ~ uniform(0, 1);

  // Likelihood
  // Forward algorithm derived from Stan Modeling Language
  // User's Guide  and Reference Manual
  for (i in 1:M) {
    // Make sure that all M individuals are in state 1 at t=1
    gam[1, 1] <- 1.0;
    gam[1, 2] <- 0.0;
    gam[1, 3] <- 0.0;

    for (t in 2:n_occasions) {
      for (k in 1:3) {
        for (j in 1:3)
          acc[j] <- gam[t - 1, j] * ps[j, i, t - 1, k]
            * po[k, i, t - 1, y[i, t]];
        gam[t, k] <- sum(acc);
      }
    }
    increment_log_prob(log(sum(gam[n_occasions])));
  }
}

generated quantities {
  real<lower=0,upper=1> psi;              // Inclusion probability
  real<lower=0,upper=1> b[n_occasions-1]; // Entry probability
  int<lower=0> Nsuper;                    // Superpopulation size
  int<lower=0> N[n_occasions-1];          // Actual population size
  int<lower=0> B[n_occasions-1];          // Number of entries
  int<lower=1,upper=3> z[M, n_occasions]; // Latent state

  // Generate z[]
  for (i in 1:M) {
    z[i, 1] <- 1;
    for (t in 2:n_occasions)
      z[i, t] <- categorical_rng(ps[z[i, t - 1], i, t - 1]);
  }

  // Calculate derived population parameters
  {
    real qgamma[n_occasions-1];
    real cprob[n_occasions-1];
    int al[M, n_occasions-1];
    int d[M, n_occasions-1];
    int alive[M];
    int w[M];

    for (t in 1:(n_occasions - 1)) {
      qgamma[t] <- 1.0 - gamma[t];
    }
    cprob[1] <- gamma[1];
    for (t in 2:(n_occasions - 1)) {
      cprob[t] <- gamma[t] * prod(qgamma[1:(t - 1)]);
    } //t
    psi <- sum(cprob);
    for (t in 1:(n_occasions - 1)) {
      b[t] <- cprob[t] / psi;
    } //t

    for (i in 1:M) {
      for (t in 2:n_occasions)
        al[i, t - 1] <- (z[i, t] == 2);
      for (t in 1:(n_occasions - 1)) {
        if (z[i, t] - al[i, t] == 0)
          d[i, t] <- 1;
        else
          d[i, t] <- 0;
      } //t
      alive[i] <- sum(al[i]);
    } //i

    for (t in 1:(n_occasions - 1)) {
      N[t] <- sum(al[1:M, t]);
      B[t] <- sum(d[1:M, t]);
    } //t
    for (i in 1:M)
      w[i] <- (alive[i] != 0);
    Nsuper <- sum(w);
  }
}
