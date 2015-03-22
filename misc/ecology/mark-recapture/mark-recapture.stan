data {
  int<lower=0> M;                  // marked
  int<lower=0> C;                  // captured
  int<lower=0,upper=min(M,C)> R;   // recaptured
}
parameters {
  real<lower=(C - R + M)> N;       // population
}
model {
  // N has improper default prior
  R ~ binomial(C, M / N);
}
