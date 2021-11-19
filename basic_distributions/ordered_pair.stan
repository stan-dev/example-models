transformed data {
  real L = -5;
  real H = 5;
}
parameters {
  real<lower=L, upper=H> a;
  real<lower=a, upper=H> b;
}
model {
  //    a ~ uniform(L, b);
  //    b ~ uniform(a, H);
}
