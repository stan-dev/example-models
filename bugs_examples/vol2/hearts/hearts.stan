/**
 * Hearts: a mixture model for count data
 * http://www.openbugs.net/Examples/Hearts.html
 * 
 * integrate out the binary parameters in hearts.stan.0 
 */
data {
  int<lower=0> N;
  array[N] int<lower=0> x;
  array[N] int<lower=0> y;
  array[N] int<lower=0> t;
}
parameters {
  real alpha;
  real delta;
}
transformed parameters {
  real<lower=0, upper=1> theta;
  theta = inv_logit(delta);
}
model {
  real p;
  real log1m_theta;
  p = inv_logit(alpha);
  log1m_theta = log1m(theta);
  
  alpha ~ normal(0, 100);
  delta ~ normal(0, 100);
  for (i in 1 : N) {
    // P(y_i = 0 | t_i) 
    //   = theta + (1 - theta) * (1 - p)^{t_i}
    // p(y_i | t_i) 
    //   = (1 - theta) * Binomial(y_i|t_i, p) for y_i = 1,2,...,t_i 
    
    if (y[i] == 0) {
      target += log(theta + (1 - theta) * pow(1 - p, t[i]));
    } else {
      target += log1m_theta + binomial_lpmf(y[i] | t[i], p);
    }
  }
}
generated quantities {
  real beta;
  beta = exp(alpha);
}
