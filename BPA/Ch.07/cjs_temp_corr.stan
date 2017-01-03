// This models is derived from section 12.3 of "Stan Modeling Language
// User's Guide and Reference Manual"

functions {
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k])
        return k;
    return 0;
  }

  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      // Compoud declaration was enabled in Stan 2.13
      int k = size(y_i) - k_rev;
      //      int k;
      //      k = size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }

  matrix prob_uncaptured(int nind, int n_occasions,
                         matrix p, matrix phi) {
    matrix[nind, n_occasions] chi;

    for (i in 1:nind) {
      chi[i, n_occasions] = 1.0;
      for (t in 1:(n_occasions - 1)) {
        // Compoud declaration was enabled in Stan 2.13
        int t_curr = n_occasions - t;
        int t_next = t_curr + 1;
        /*
        int t_curr;
        int t_next;

        t_curr = n_occasions - t;
        t_next = t_curr + 1;
        */
        chi[i, t_curr] = (1 - phi[i, t_curr])
                        + phi[i, t_curr]  * (1 - p[i, t_next - 1]) * chi[i, t_next];
      }
    }
    return chi;
  }
}

data {
  int<lower=0> nind;            // Number of individuals
  int<lower=2> n_occasions;     // Number of capture occasions
  int<lower=0,upper=1> y[nind, n_occasions];    // Capture-history
  int<lower=1> g;               // Number of groups
  int<lower=1,upper=g> group[nind];     // Groups
  real df;                      // Degree of freedom
  matrix[g,g] R;                // Scale matrix
}

transformed data {
  // Compoud declaration is enabled in Stan 2.13
  int n_occ_minus_1 = n_occasions - 1;
  //  int n_occ_minus_1;
  int<lower=0,upper=n_occasions> first[nind];
  int<lower=0,upper=n_occasions> last[nind];

  //  n_occ_minus_1 = n_occasions - 1;
  for (i in 1:nind)
    first[i] = first_capture(y[i]);
  for (i in 1:nind)
    last[i] = last_capture(y[i]);
}

parameters {
  vector<lower=0,upper=1>[g] mean_phi;   // Mean group-spec. survival
  real<lower=0,upper=1> p_g[g];          // Group-spec. recapture
  matrix[n_occ_minus_1, g] eta_phi;
  cov_matrix[g] Sigma;
}

transformed parameters {
  matrix<lower=0,upper=1>[nind, n_occ_minus_1] phi;
  matrix<lower=0,upper=1>[nind, n_occ_minus_1] p;
  matrix<lower=0,upper=1>[nind, n_occasions] chi;
  vector[g] mu_phi;

  for (u in 1:g)
    mu_phi[u] = logit(mean_phi[u]);

  // Constraints
  for (i in 1:nind) {
    for (t in 1:(first[i] - 1)) {
      phi[i, t] = 0;
      p[i, t] = 0;
    }
    for (t in first[i]:n_occ_minus_1) {
      phi[i, t] = inv_logit(eta_phi[t, group[i]]);
      p[i, t] = p_g[group[i]];
    }
  }

  chi = prob_uncaptured(nind, n_occasions, p, phi);
}

model {
  // Priors
  // for survival parameters
  Sigma ~ inv_wishart(df, R);

  for (t in 1:n_occ_minus_1)
    eta_phi[t] ~ multi_normal(mu_phi, Sigma);
  // Uniform priors are implicitly defined.
  //  mean_phi ~ uniform(0, 1);

  // for recapture parameters
  //  p_g ~ uniform(0, 1);

  // Likelihood
  for (i in 1:nind) {
    if (first[i] > 0) {
      for (t in (first[i] + 1):last[i]) {
        1 ~ bernoulli(phi[i, t - 1]);
        y[i, t] ~ bernoulli(p[i, t - 1]);
      }
      1 ~ bernoulli(chi[i, last[i]]);
    }
  }
}
