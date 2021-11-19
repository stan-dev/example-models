data {
  int<lower=1> I; // # of items
  int<lower=1> J; // # of respondents	
  int<lower=1> K; // # of attributes
  int<lower=1> C; // # of attribute profiles (latent classes)	
  matrix[J, I] y; // response matrix
  matrix[C, K] alpha; // attribute profile matrix
  matrix[I, C] xi; // the global attribute mastery indicator (product of alpha^q-element)
}
parameters {
  simplex[C] nu; // probabilities of latent class membership
  array[I] real<lower=0, upper=1> slip; // slip parameter	
  array[I] real<lower=0, upper=1> guess; // guess parameter
}
transformed parameters {
  vector[C] log_nu = log(nu);
}
model {
  array[C] real ps; // temp for log component densities
  matrix[I, C] pi;
  array[I] real log_items;
  slip ~ beta(5, 25);
  guess ~ beta(5, 25);
  for (c in 1 : C) {
    for (i in 1 : I) {
      pi[i, c] = (1 - slip[i]) ^ xi[i, c] * guess[i] ^ (1 - xi[i, c]);
    }
  }
  for (j in 1 : J) {
    for (c in 1 : C) {
      for (i in 1 : I) {
        log_items[i] = y[j, i] * log(pi[i, c])
                       + (1 - y[j, i]) * log(1 - pi[i, c]);
      }
      ps[c] = log_nu[c] + sum(log_items);
    }
    target += log_sum_exp(ps);
  }
}
generated quantities {
  matrix[J, C] prob_resp_class; // posterior probabilities of respondent j being in latent class c 
  matrix[J, K] prob_resp_attr; // posterior probabilities of respondent j being a master of attribute k 
  matrix[I, C] pi;
  array[I] real log_items;
  row_vector[C] prob_joint;
  array[C] real prob_attr_class;
  for (c in 1 : C) {
    for (i in 1 : I) {
      pi[i, c] = (1 - slip[i]) ^ xi[i, c] * guess[i] ^ (1 - xi[i, c]);
    }
  }
  for (j in 1 : J) {
    for (c in 1 : C) {
      for (i in 1 : I) {
        log_items[i] = y[j, i] * log(pi[i, c])
                       + (1 - y[j, i]) * log(1 - pi[i, c]);
      }
      prob_joint[c] = nu[c] * exp(sum(log_items));
    }
    prob_resp_class[j] = prob_joint / sum(prob_joint);
  }
  for (j in 1 : J) {
    for (k in 1 : K) {
      for (c in 1 : C) {
        prob_attr_class[c] = prob_resp_class[j, c] * alpha[c, k];
      }
      prob_resp_attr[j, k] = sum(prob_attr_class);
    }
  }
}
