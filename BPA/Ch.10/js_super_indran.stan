// JS model using the superpopulation parameterization

functions {
  // These functions are derived from Section 12.3 of
  // Stan Modeling Language User's Guide and Reference Manual

  /**
   * Return a integer value of first capture occasion
   *
   * @param y_i Integer array of capture history
   *
   * @return Integer value of first capture occasion
   */
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k])
        return k;
    return 0;
  }

  /**
   * Return a integer value of last capture occasion
   *
   * @param y_i Integer array of capture history
   *
   * @return Integer value of last capture occasion
   */
  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      int k = size(y_i) - k_rev;

      if (y_i[k])
        return k;
    }
    return 0;
  }

  /**
   * Return a matrix of uncaptured probabilities
   *
   * @param p   Matrix of detection probabilities for each individual
   *            and capture occasion
   * @param phi Matrix of survival probabilities for each individual
   *            and capture occasion
   *
   * @return Matrix of uncaptured probabilities
   */
  matrix prob_uncaptured(matrix p, matrix phi) {
    int n_ind = rows(p);
    int n_occasions = cols(p);
    matrix[n_ind, n_occasions] chi;

    for (i in 1:n_ind) {
      chi[i, n_occasions] = 1.0;
      for (t in 1:(n_occasions - 1)) {
        int t_curr = n_occasions - t;
        int t_next = t_curr + 1;

        chi[i, t_curr] = (1 - phi[i, t_curr])
          + phi[i, t_curr] * (1 - p[i, t_next]) * chi[i, t_next];
      }
    }
    return chi;
  }

  /**
   * Calculate log likelihood of a Jolly-Seber model
   * under the superpopulation parameterization
   *
   * @param y     Integer array of capture history
   * @param first Integer array of first capture occasions
   * @param last  Integer array of last capture occasions
   * @param p     Matrix of detection probabilities
   * @param phi   Matrix of survival probabilities
   * @param psi   Real value of inclusion probability
   * @param nu    Vector of entry probabilities
   * @param chi   Matrix of uncapture probabilities
   */
  void js_super_lp(int[,] y, int[] first, int[] last,
                   matrix p, matrix phi, real psi,
                   vector nu, matrix chi) {
    int n_ind = dims(y)[1];
    int n_occasions = dims(y)[2];
    vector[n_occasions] qnu = 1.0 - nu;

    for (i in 1:n_ind) {
      vector[n_occasions] qp = 1.0 - p[i]';

      if (first[i]) { // Observed
        // Included
        1 ~ bernoulli(psi);

        // Until first capture
        if (first[i] == 1) {
          1 ~ bernoulli(nu[1] * p[i, 1]);
        } else {  // first[i] >= 2
          vector[first[i]] lp;

          // Entered at 1st occasion
          lp[1] = bernoulli_lpmf(1 | nu[1])
                 + bernoulli_lpmf(1 | prod(qp[1:(first[i] - 1)]))
                 + bernoulli_lpmf(1 | prod(phi[i, 1:(first[i] - 1)]))
                 + bernoulli_lpmf(1 | p[i, first[i]]);
          // Entered at t-th occasion (1 < t < first[i])
          for (t in 2:(first[i] - 1))
            lp[t] = bernoulli_lpmf(1 | prod(qnu[1:(t - 1)]))
                   + bernoulli_lpmf(1 | nu[t])
                   + bernoulli_lpmf(1 | prod(qp[t:(first[i] - 1)]))
                   + bernoulli_lpmf(1 | prod(phi[i, t:(first[i] - 1)]))
                   + bernoulli_lpmf(1 | p[i, first[i]]);
          lp[first[i]] = bernoulli_lpmf(1 | prod(qnu[1:(first[i] - 1)]))
                        + bernoulli_lpmf(1 | nu[first[i]])
                        + bernoulli_lpmf(1 | p[i, first[i]]);
          target += log_sum_exp(lp);
        }
        // Until last capture
        for (t in (first[i] + 1):last[i]) {
          1 ~ bernoulli(phi[i, t - 1]);   // Survived
          y[i, t] ~ bernoulli(p[i, t]);   // Capture/Non-capture
        }
        // Subsequent occasions
        1 ~ bernoulli(chi[i, last[i]]);
      } else {          // Never observed
        vector[n_occasions+1] lp;

        // Entered at 1st occasion, but never captured
        lp[1] = bernoulli_lpmf(1 | psi)
               + bernoulli_lpmf(1 | nu[1])
               + bernoulli_lpmf(0 | p[i, 1])
               + bernoulli_lpmf(1 | chi[i, 1]);
        // Entered at t-th occation (t > 1), but never captured
        for (t in 2:n_occasions)
          lp[t] = bernoulli_lpmf(1 | psi)
                 + bernoulli_lpmf(1 | prod(qnu[1:(t - 1)]))
                 + bernoulli_lpmf(1 | nu[t])
                 + bernoulli_lpmf(0 | p[i, t])
                 + bernoulli_lpmf(1 | chi[i, t]);
        // Never captured
        lp[n_occasions + 1] = bernoulli_lpmf(0 | psi);
        target += log_sum_exp(lp);
      }
    }
  }
}

data {
  int<lower=0> M;                              // Augmented sample size
  int<lower=0> n_occasions;                    // Number of capture occasions
  int<lower=0,upper=1> y[M, n_occasions];      // Augmented capture-history
}

transformed data {
  int<lower=0,upper=n_occasions> first[M];
  int<lower=0,upper=n_occasions> last[M];

  for (i in 1:M)
    first[i] = first_capture(y[i]);
  for (i in 1:M)
    last[i] = last_capture(y[i]);
}

parameters {
  real<lower=0,upper=1> mean_phi;             // Mean survival
  real<lower=0,upper=1> mean_p;               // Mean capture
  real<lower=0,upper=1> psi;                  // Inclusion probability
  vector<lower=0>[n_occasions] beta;
  vector[M] epsilon;
  real<lower=0,upper=5> sigma;
  // In case of using a weakly informative prior
  //  real<lower=0> sigma;
}

transformed parameters {
  matrix<lower=0,upper=1>[M, n_occasions - 1] phi;
  matrix<lower=0,upper=1>[M, n_occasions] p;
  simplex[n_occasions] b;                     // Entry probability
  vector<lower=0,upper=1>[n_occasions] nu;
  matrix<lower=0,upper=1>[M, n_occasions] chi;

  // Constraints
  phi = rep_matrix(mean_phi, M, n_occasions - 1);
  for (t in 1:n_occasions)
    p[, t] = inv_logit(logit(mean_p) + epsilon);


  // Dirichlet prior for entry probabilities
  // beta ~ gamma(1, 1);  // => model block
  b = beta / sum(beta);

  // Convert entry probs to conditional entry probs
  {
    real cum_b = b[1];

    nu[1] = b[1];
    for (t in 2:(n_occasions - 1)) {
      nu[t] = b[t] / (1.0 - cum_b);
      cum_b = cum_b + b[t];
    }
    nu[n_occasions] = 1.0;
  }

  // Uncaptured probability
  chi = prob_uncaptured(p, phi);
}

model {
  // Priors
  // Uniform priors are implicitly defined.
  // In case of using a weakly informative prior on sigma
  //  sigma ~ normal(2.5, 1.25);
  epsilon ~ normal(0, sigma);
  beta ~ gamma(1, 1);

  // Likelihood
  js_super_lp(y, first, last, p, phi, psi, nu, chi);
}

generated quantities {
  real<lower=0> sigma2;
  int<lower=0> Nsuper;                    // Superpopulation size
  int<lower=0> N[n_occasions];            // Actual population size
  int<lower=0> B[n_occasions];            // Number of entries
  int<lower=0,upper=1> z[M, n_occasions]; // Deflated latent state

  sigma2 = square(sigma);

  // Generate w[] and z[]
  for (i in 1:M) {
    int q = 1;
    if (bernoulli_rng(psi)) {      // Included
      z[i, 1] = bernoulli_rng(nu[1]);
      for (t in 2:n_occasions) {
        q = q * (1 - z[i, t - 1]);
        z[i, t] = bernoulli_rng(z[i, t - 1] * phi[i, t - 1]
                                 + q * nu[t]);
      }
    } else {                       // Not included
      z[i, ] = rep_array(0, n_occasions);
    }
  }

  // Calculate derived population parameters
  {
    int recruit[M, n_occasions] = rep_array(0, M, n_occasions);
    int Nind[M];
    int Nalive[M];

    for (i in 1:M) {
      int f = first_capture(z[i, ]);

      if (f > 0)
        recruit[i, f] = 1;
    }
    for (t in 1:n_occasions) {
      N[t] = sum(z[, t]);
      B[t] = sum(recruit[, t]);
    }
    for (i in 1:M) {
      Nind[i] = sum(z[i]);
      Nalive[i] = (Nind[i] > 0);
    }
    Nsuper = sum(Nalive);
  }
}
