data {
  int<lower=1> K; // num categories
  int<lower=1> V; // num words
  int<lower=0> T; // num supervised items
  int<lower=1> T_unsup; // num unsupervised items
  array[T] int<lower=1, upper=V> w; // words
  array[T] int<lower=1, upper=K> z; // categories
  array[T_unsup] int<lower=1, upper=V> u; // unsup words
  vector<lower=0>[K] alpha; // transit prior
  vector<lower=0>[V] beta; // emit prior
}
parameters {
  array[K] simplex[K] theta; // transit probs
  array[K] simplex[V] phi; // emit probs
}
model {
  for (k in 1 : K) {
    theta[k] ~ dirichlet(alpha);
  }
  for (k in 1 : K) {
    phi[k] ~ dirichlet(beta);
  }
  for (t in 1 : T) {
    w[t] ~ categorical(phi[z[t]]);
  }
  for (t in 2 : T) {
    z[t] ~ categorical(theta[z[t - 1]]);
  }
  
  {
    // forward algorithm computes log p(u|...)
    array[K] real acc;
    array[T_unsup, K] real gamma;
    for (k in 1 : K) {
      gamma[1, k] = log(phi[k, u[1]]);
    }
    for (t in 2 : T_unsup) {
      for (k in 1 : K) {
        for (j in 1 : K) {
          acc[j] = gamma[t - 1, j] + log(theta[j, k]) + log(phi[k, u[t]]);
        }
        gamma[t, k] = log_sum_exp(acc);
      }
    }
    target += log_sum_exp(gamma[T_unsup]);
  }
}
generated quantities {
  array[T_unsup] int<lower=1, upper=K> y_star;
  real log_p_y_star;
  {
    // Viterbi algorithm
    array[T_unsup, K] int back_ptr;
    array[T_unsup, K] real best_logp;
    for (k in 1 : K) {
      best_logp[1, K] = log(phi[k, u[1]]);
    }
    for (t in 2 : T_unsup) {
      for (k in 1 : K) {
        best_logp[t, k] = negative_infinity();
        for (j in 1 : K) {
          real logp;
          logp = best_logp[t - 1, j] + log(theta[j, k]) + log(phi[k, u[t]]);
          if (logp > best_logp[t, k]) {
            back_ptr[t, k] = j;
            best_logp[t, k] = logp;
          }
        }
      }
    }
    log_p_y_star = max(best_logp[T_unsup]);
    for (k in 1 : K) {
      if (best_logp[T_unsup, k] == log_p_y_star) {
        y_star[T_unsup] = k;
      }
    }
    for (t in 1 : (T_unsup - 1)) {
      y_star[T_unsup - t] = back_ptr[T_unsup - t + 1, y_star[T_unsup - t + 1]];
    }
  }
}
