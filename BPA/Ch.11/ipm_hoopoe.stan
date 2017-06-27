// Integrated population model

functions {
  /**
   * Return log probability of Poisson distribution.
   * outcome n may be a real value; for compatibility with Win/OpenBUGS.
   *
   * @param n      Outcome
   * @param lambda Mean
   *
   * @return Log probability
   */
  real real_poisson_lpdf(real n, real lambda) {
    real lp;

    if (lambda < 0) {
      reject("lambda must be non-negative; found lambda=", lambda);
    } else if (n < 0) {
      reject("n must not be negative; found n=", n);
    } else {
      lp = n * log(lambda) - lambda - lgamma(n + 1);
    }
    return lp;
  }

  /**
   * Return log probability of binomial distribution.
   * outcome n may be a real value; for compatibility with Win/OpenBUGS.
   *
   * @param n     Outcome
   * @param N     Size
   * @param theta Probability
   *
   * @return Log probability
   */
  real real_binomial_lpdf(real n, real N, real theta) {
    real lp;

    if (N < 0) {
      reject("N must be non-negative; found N=", N);
    } else if (theta < 0 || theta > 1) {
      reject("theta must be in [0,1]; found theta=", theta);
    } else if (n < 0 || n > N) {
      reject("n must be in [0,N]; found n=", n);
    } else {
      lp = lchoose(N, n) + n * log(theta) + (N - n) * log(1 - theta);
    }
    return lp;
  }

  /**
   * Return m-array  cell probabilities for juveniles
   *
   * @param nyears Number of years
   * @param phij   Survival probability
   * @param p      Recapture probability
   *
   * @return m-array cell probabilities for juveniles
   */
  vector[] marray_juveniles(int nyears, vector phij, vector phia, vector p) {
    int ny_minus_1 = nyears - 1;
    vector[nyears] pr_j[ny_minus_1];
    vector[ny_minus_1] q = 1 - p;
    real prod_phi;
    real prod_q;

    for (t in 1:ny_minus_1) {
      // Main diagonal
      pr_j[t, t] = phij[t] * p[t];

      // Above main diagonal
      prod_phi = 1;
      prod_q = 1;
      for (j in (t + 1):ny_minus_1) {
        prod_phi = prod_phi * phia[j];
        prod_q = prod_q * q[j - 1];
        pr_j[t, j] = phij[t] * prod_phi * prod_q * p[j];
      }
      // Below main diagonal
      for (j in 1:(t - 1))
        pr_j[t, j] = 0;

      // Last column: probability of non-recapture
      pr_j[t, nyears] = 1 - sum(pr_j[t, 1:(nyears - 1)]);
    }
    return pr_j;
  }

  /**
   * Return m-array cell probabilities for adults
   *
   * @param nyears Number of years
   * @param phia   Survival probability
   * @param p      Recapture probability
   *
   * @return m-array cell probabilities for adults
   */
  vector[] marray_adults(int nyears, vector phia, vector p) {
    int ny_minus_1 = nyears - 1;
    vector[nyears] pr_a[ny_minus_1];
    vector[nyears-1] q = 1 - p;
    real prod_phi;
    real prod_q;

    for (t in 1:(nyears - 1)) {
      // Main diagonal
      pr_a[t, t] = phia[t] * p[t];

      // Above main diagonal
      prod_phi = phia[t];
      prod_q = 1;
      for (j in (t + 1):(nyears - 1)) {
        prod_phi = prod_phi * phia[j];
        prod_q = prod_q * q[j - 1];
        pr_a[t, j] = prod_phi * prod_q * p[j];
      }

      // Below main diagonal
      for (j in 1:(t - 1))
        pr_a[t, j] = 0;

      // Last column
      pr_a[t, nyears] = 1 - sum(pr_a[t, 1:(nyears - 1)]);
    }
    return pr_a;
  }
}

data {
  int nyears;                        // Number of years
  int y[nyears];                     // Population counts
  int J[nyears - 1];                 // Productivity data
  int R[nyears - 1];
  int marray_j[nyears - 1, nyears];  // Data in m-array format for juveniles
  int marray_a[nyears - 1, nyears];  // Data in m-array format for adults
}

transformed data {
  int ny_minus_1 = nyears - 1;
}

parameters {
  vector<lower=0>[nyears] N1;      // Number of 1-year old individuals
  vector<lower=0>[nyears] NadSurv; // Number of Adults >= 2 years
  vector<lower=0>[nyears] Nadimm;  // Number of immigrants
  real l_mphij;
  real l_mphia;
  real l_mfec;
  real l_mim;
  real l_p;
  vector[ny_minus_1] epsilon_phij_raw;
  vector[ny_minus_1] epsilon_phia_raw;
  vector[ny_minus_1] epsilon_fec_raw;
  vector[ny_minus_1] epsilon_im_raw;
  real<lower=0> sig_phij;
  real<lower=0> sig_phia;
  real<lower=0> sig_fec;
  real<lower=0> sig_im;
}

transformed parameters {
  vector[ny_minus_1] epsilon_phij;
  vector[ny_minus_1] epsilon_phia;
  vector[ny_minus_1] epsilon_fec;
  vector[ny_minus_1] epsilon_im;
  vector<lower=0,upper=1>[ny_minus_1] phij;  // Juv. apparent survival
  vector<lower=0,upper=1>[ny_minus_1] phia;  // Adult apparent survival
  vector<lower=0>[ny_minus_1]         f;     // Productivity
  vector<lower=0>[ny_minus_1]         omega; // Immigration
  vector<lower=0,upper=1>[ny_minus_1] p;     // Recapture probability
  vector<lower=0>[nyears] Ntot;              // Number of total individuals
  vector<lower=0>[ny_minus_1] rho;
  simplex[nyears] pr_j[ny_minus_1];
  simplex[nyears] pr_a[ny_minus_1];

  // Distribution of error terms
  epsilon_phij = sig_phij * epsilon_phij_raw;
  epsilon_phia = sig_phia * epsilon_phia_raw;
  epsilon_fec = sig_fec * epsilon_fec_raw;
  epsilon_im = sig_im * epsilon_im_raw;

  // Constrain parameters
  for (t in 1:(nyears - 1)) {
    phij[t] = inv_logit(l_mphij + epsilon_phij[t]);
    phia[t] = inv_logit(l_mphia + epsilon_phia[t]);
    f[t] = exp(l_mfec + epsilon_fec[t]);
    omega[t] = exp(l_mim + epsilon_im[t]);
    p[t] = inv_logit(l_p);
  }

  // Total number of individuals
  Ntot = NadSurv + Nadimm + N1;

  // m-array
  pr_j = marray_juveniles(nyears, phij, phia, p);
  pr_a = marray_adults(nyears, phia, p);

  // Productivity
  for (t in 1:ny_minus_1)
    rho[t] = R[t] * f[t];
}

model {
  // Priors
  // Initial population sizes
  // Constraints ensure truncated normal (> 0)
  N1[1] ~ normal(100, 100);
  NadSurv[1] ~ normal(100, 100);
  Nadimm[1] ~ normal(100, 100);

  // Mean demographic parameters (on appropriate scale)
  l_mphij ~ normal(0, 100);
  l_mphia ~ normal(0, 100);
  l_mfec ~ normal(0, 100);
  l_mim ~ normal(0, 100);
  l_p ~ normal(0, 100);

  // Improper flat priors are implicitly used on
  // scale parameters, sig_phij, sig_phia, sig_fec and sig_im.

  // Efficient form of priors
  epsilon_phij_raw ~ normal(0, 1);
  epsilon_phia_raw ~ normal(0, 1);
  epsilon_fec_raw ~ normal(0, 1);
  epsilon_im_raw ~ normal(0, 1);

  // Likelihood for population population count data (state-space model)
  // System process
  for (t in 2:nyears) {
    real mean1 = 0.5 * f[t - 1] * phij[t - 1] * Ntot[t - 1];
    real mpo = Ntot[t - 1] * omega[t - 1];

    N1[t] ~ real_poisson(mean1);
    NadSurv[t] ~ real_binomial(Ntot[t - 1], phia[t - 1]);
    Nadimm[t] ~ real_poisson(mpo);
  }

  // Observation process
  y ~ poisson(Ntot);

  // Likelihood for capture-recapture data: CJS model (2 age classes)
  // Multinomial likelihood
  for (t in 1:(nyears - 1)) {
    marray_j[t] ~ multinomial(pr_j[t]);
    marray_a[t] ~ multinomial(pr_a[t]);
  }

  // Likelihood for productivity data: Poisson regression
  J ~ poisson(rho);
}

generated quantities {
  real<lower=0,upper=1> mphij = inv_logit(l_mphij); // Mean juv. survival prob.
  real<lower=0,upper=1> mphia = inv_logit(l_mphia); // Mean adult survival prob.
  real<lower=0> mfec = exp(l_mfec);    // Mean productivity
  real<lower=0> mim = exp(l_mim);      // Mean immigration rate
  vector<lower=0>[ny_minus_1] lambda;  // Population growth rate
  vector[ny_minus_1] logla;
  real<lower=0> mlam;

  // Population growth rate
  lambda[1:(nyears - 1)] = Ntot[2:nyears] ./ Ntot[1:(nyears - 1)];
  logla = log(lambda);
  // Geometric mean
  mlam = exp(1.0 / ny_minus_1 * sum(logla));
}
