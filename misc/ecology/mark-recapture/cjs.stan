/*
 * Cormack-Jolly-Server (CJS) model.
 *
 * Key is that log probability is conditioned on first observation for
 * an individual.  See Matt Schofield's Ph.D. thesis for a nice
 * overview of CJS and generalizations:
 *
 * Schofield, Matthew R.  2007.  Hierarchical Capture-Recapture
 * Models.  Dept. of Statistics, University of Otago, Dunedin.
 *
 */
data {
  int<lower=0> history[7];
}
parameters {
  real<lower=0,upper=1> phi[2]; // phi[t] survival from t to t+1
  real<lower=0,upper=1> p[3];   // p[t]   probability of capture at t+1
}
transformed parameters {
  //  chi[t] probability a sighting at t is last sighting
  real<lower=0,upper=1> chi[3];  
  chi[3] <- 1;
  chi[2] <- (1 - phi[2]) + phi[2] * (1 - p[3]);
  chi[1] <- (1 - phi[1]) + phi[1] * (1 - p[2]) * chi[2];
}
model {
  increment_log_prob(history[7] * (log(phi[1]) + log(p[2]) + log(phi[2]) + log(p[3])));
  increment_log_prob(history[6] * (log(phi[1]) + log(p[2]) + log(chi[2])));
  increment_log_prob(history[5] * (log(phi[1]) + log1m(p[2]) + log(phi[2]) + log(p[3])));
  increment_log_prob(history[4] * (log(chi[1])));
  increment_log_prob(history[3] * (log(phi[2]) + log(p[3])));
  increment_log_prob(history[2] * log(chi[2]));
  // history[1] provides no information
}
generated quantities {
  real<lower=0,upper=1> beta3;
  beta3 <- phi[2] * p[3];
}
