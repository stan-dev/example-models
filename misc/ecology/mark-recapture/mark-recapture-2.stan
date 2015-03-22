data {
  int<lower=0> M;                       // marked
  int<lower=0> C;                       // captured
  int<lower=0,upper=min(M,C)> R;        // recaptured
}
transformed data {
  real theta_max;
  theta_max <- M;         
  theta_max <- theta_max / (C - R + M);
}
parameters {
  real<lower=0,upper=theta_max> theta;  // proportion marked
}
model {
  R ~ binomial(C, theta);
}
generated quantities {
  real<lower=(C - R + M)> N;            // population
  N <- M / theta;
}

