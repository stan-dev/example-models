/*
 * mark-recpapture-2 prior:
 *
 *   theta ~ uniform(0, M / (C - R + M))
 * 
 * reparameterized in temrs of N:
 *
 *   (M / N) ~ uniform(0, M / (C - R + M));
 * 
 * requires Jacobian adjustment:
 *
 *   function:      f(N) -> M / N
 *   derivative:    f'(N)  =  -M / N^2
 *   Jacobian adj:  abs(f'(N)) = M / N^2
 *   log Jac adj:   log(abs(f'(N)) = log(M) - 2 * log(N)
 *   drop const:    -2 * log(N)
 */

data {
  int<lower=0> M;
  int<lower=0> C;
  int<lower=0,upper=min(M,C)> R;
}
transformed data {
  real theta_max;
  theta_max <- M;         
  theta_max <- theta_max / (C - R + M);
}
parameters {
  real<lower=(C - R + M)> N;
}
transformed parameters {
  real<lower=0,upper=theta_max> theta;
  theta <- M / N;
}
model {
  increment_log_prob(-2 * log(N));
  R ~ binomial(C, theta);
}
