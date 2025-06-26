// #### Notes to Stan model #######################################################
// ## Stan has hypergeometric distribution implemented, so in one way the code
// ## is more intuitive. On the other hand, Stan can't sample discrete parameters,
// ## therefore we have to increment log probability manually (as we did in 
// ## Survey example).  
// ################################################################################

// Planes
data { 
  int<lower=1> x;
  int<lower=1> n;
  int<lower=0, upper=n> k;
  int<lower=x> tmax;
}

transformed data {
  int<lower=x> tmin;
    
  tmin = x + n - k;
}

transformed parameters {
  vector[tmax] lp_parts;

  for (t in 1:tmax) 
    if (t < tmin)
      lp_parts[t] = log(1.0 / tmax) + negative_infinity();  // Zero probability
    else
      lp_parts[t] = log(1.0 / tmax) + hypergeometric_lpmf(k | n, x, t - x);
}

model {
  target += log_sum_exp(lp_parts);
}

generated quantities {
  int<lower=tmin, upper=tmax> t;
  simplex[tmax] tp;

  tp = softmax(lp_parts);
  t = categorical_rng(tp);
}
