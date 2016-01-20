functions {
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k])
        return k;
    return 0;
  }

  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      int k;
      k <- size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }

  real cell_prob(int n_occasions, row_vector s, row_vector r,
                   int first, int last) {
    vector[n_occasions] pr;

    pr[first] <- (1 - s[first]) * r[first];

    for (j in (first + 1):n_occasions - 1)
      pr[j] <- prod(segment(s, first, j - first)) * (1 - s[j]) * r[j];

    for (j in 1:(first - 1))
      pr[j] <- 0;

    for (t in 1:n_occasions - 1)
      pr[n_occasions] <- 1 - sum(head(pr, n_occasions - 1));

    return pr[last];
  }
}

data {
  int<lower=0> nind;
  int<lower=0> n_occasions;
  int<lower=0,upper=1> y[nind, n_occasions];
}

transformed data {
  int<lower=0,upper=n_occasions> first[nind];
  int<lower=0,upper=n_occasions> last[nind];

  for (i in 1:nind)
    first[i] <- first_capture(y[i]);
  for (i in 1:nind)
    last[i] <- last_capture(y[i]);
}

parameters {
  real<lower=0,upper=1> mean_s;     // Mean survival
  real<lower=0,upper=1> mean_r;     // Mean recapture
}

transformed parameters {
  matrix<lower=0,upper=1>[nind, n_occasions - 1] s;
  matrix<lower=0,upper=1>[nind, n_occasions - 1] r;

  // Constraints
  for (i in 1:nind) {
    for (t in 1:(first[i] - 1)) {
      s[i, t] <- 0;
      r[i, t] <- 0;
    }
    for (t in first[i]:n_occasions - 1) {
      s[i, t] <- mean_s;
      r[i, t] <- mean_r;
    }
  }
}

model {
  // Priors
  mean_s ~ uniform(0, 1);
  mean_r ~ uniform(0, 1);

  // Likelihood 
  for (i in 1:nind){
    real pr;

    if (first[i] > 0) {
      if (last[i] > first[i]) // Recovered
        pr <- cell_prob(n_occasions, s[i], r[i], first[i], last[i] - 1);
      else                    // Never Recovered
        pr <- cell_prob(n_occasions, s[i], r[i], first[i], n_occasions);
      1 ~ bernoulli(pr);
    }
  }
}
