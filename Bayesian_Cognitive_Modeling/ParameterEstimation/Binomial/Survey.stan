/* 
 * Notes to Stan model 
 * --------------------
 * This model code is more difficult to understand in Stan implementation since 
 * Stan is unable to sample discrete parameters. This may change in the future. 
 * For better understanding read Stan manual chapter "Mixture Modeling" (p.68 in 
 * version 2.4.0) first.
 */

// Inferring Return Rate and Number of Surveys from Observed Returns
data { 
  int<lower=0> nmax;
  int<lower=0> m;
  int<lower=0,upper=nmax> k[m];
}
transformed data {
  int<lower=0> nmin;  // Minimal possible n
  
  nmin <- max(k);
}
parameters {
  real<lower=0,upper=1> theta;
}
transformed parameters {
  vector[nmax] lp_parts;  // Log probability for each n

  // First part of the trick for mixture model
  for (n in 1:nmax)
    if (n < nmin)
      lp_parts[n] <- log(1.0 / nmax) + negative_infinity();  // Zero probability
    else
      lp_parts[n] <- log(1.0 / nmax) + binomial_log(k, n, theta); 
}
model {
  // Second part of the trick for mixture model
  increment_log_prob(log_sum_exp(lp_parts));
}
generated quantities {
  int<lower=1,upper=nmax> n;
  simplex[nmax] prob_n;
  
  // Transforming lp_parts to probabilities of each n
  prob_n <- softmax(lp_parts);
  n <- categorical_rng(prob_n);
}