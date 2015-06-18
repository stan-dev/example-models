data {
  int<lower=1> I;               // # items
  int<lower=1> J;               // # students
  int<lower=0,upper=1> y[I,J];  // boolean correctness of student j on item i for n

}
parameters {
  vector[J] alpha;             // ability for student j
  vector[I] beta;              // difficulty for item i
  vector<lower=0>[I] delta;    // discriminativeness for item i
}
model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 5);
  delta ~ gamma(2, 2);
  for (i in 1:I)
    y[i] ~ bernoulli_logit(delta[i] * (alpha - beta[i]));
}
