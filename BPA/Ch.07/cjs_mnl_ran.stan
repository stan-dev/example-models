data {
  int<lower=0> n_occasions;     // Number of capture occasions
  int<lower=0> marr[n_occasions - 1, n_occasions];      // m-array
}

transformed data {
  // Compoud declaration is enabled in Stan 2.13
  int n_occ_minus_1 = n_occasions - 1;
  //  int n_occ_minus_1;
  int r[n_occasions];

  //  n_occ_minus_1 = n_occasions - 1;
  for (t in 1:n_occ_minus_1)
    r[t] = sum(marr[t]);
}

parameters {
  real<lower=0,upper=1> mean_phi;       // Mean survival
  real<lower=0,upper=1> mean_p;         // Mean recapture
  real<lower=0> sigma;
  vector[n_occ_minus_1] epsilon;
}

transformed parameters {
  real mu;
  vector<lower=0,upper=1>[n_occ_minus_1] phi;
  vector<lower=0,upper=1>[n_occ_minus_1] p;
  vector<lower=0,upper=1>[n_occ_minus_1] q;
  simplex[n_occasions] pr[n_occ_minus_1];

  mu = logit(mean_phi);

  // Constraints
  // inv_logit was vectorized in Stan 2.13
  phi = inv_logit(mu + epsilon);
  //  for (t in 1:n_occ_minus_1)
  //    phi[t] = inv_logit(mu + epsilon[t]);
  p = rep_vector(mean_p, n_occ_minus_1);

  q = 1.0 - p;             // Probability of non-recapture

  // Define the cell probabilities of the m-arrays
  // Main diagonal
  for (t in 1:n_occ_minus_1) {
    pr[t, t] = phi[t] * p[t];

    // Above main diagonal
    for (j in (t + 1):n_occ_minus_1)
      pr[t, j] = prod(phi[t:j]) * prod(q[t:(j - 1)]) * p[j];

    // Below main diagonal
    pr[t, :(t - 1)] = rep_vector(0, t - 1);
  }

  // Last column: probability of non-recapture
  for (t in 1:n_occ_minus_1)
    pr[t, n_occasions] = 1 - sum(pr[t, :n_occ_minus_1]);
}

model {
  // Priors
  // Uniform priors are implicitly defined.
  //  mean_phi ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);
  //  sigma ~ uniform(0, 5);
  // In case a weakly informative prior is used
  //  sigma ~ normal(2.5, 1.25);
  epsilon ~ normal(0, sigma);

  // Define the multinomial likelihood
  for (t in 1:(n_occ_minus_1))
    marr[t] ~ multinomial(pr[t]);
}

generated quantities {
  real<lower=0> sigma2;
  real<lower=0> sigma2_real;
  vector[n_occasions] expmarr[n_occ_minus_1];
  int marr_new[n_occ_minus_1, n_occasions];
  matrix[n_occ_minus_1, n_occasions] E_org;
  matrix[n_occ_minus_1, n_occasions] E_new;
  real fit;
  real fit_new;

  sigma2 = square(sigma);

  // Temporal variance on real scale
  sigma2_real = sigma2 * square(mean_phi * (1 - mean_phi));

  // Assess model fit using Freeman-Tukey statistic
  // Compute fit statistics for observed data
  for (t in 1:n_occ_minus_1) {
    expmarr[t] = r[t] * pr[t];
    for (j in 1:n_occasions)
      E_org[t, j] = square((sqrt(marr[t, j]) - sqrt(expmarr[t][j])));
   }

  // Generate replicate data and compute fit stats from them
  for (t in 1:n_occ_minus_1) {
    marr_new[t] = multinomial_rng(pr[t], r[t]);
    for (j in 1:n_occasions)
      E_new[t, j] = square((sqrt(marr_new[t, j]) - sqrt(expmarr[t][j])));
  }
  fit = sum(E_org);
  fit_new = sum(E_new);
}
