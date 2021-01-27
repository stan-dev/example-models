parameters {
  real<lower=-1,upper=1> y;
}
model {
  target += log1m(fabs(y));
}
