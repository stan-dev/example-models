// IRT 1PL (RASCH) MODEL 

data {
  int<lower=1> J; // number of students
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of observations
  array[N] int<lower=1, upper=J> jj; // student for observation n
  array[N] int<lower=1, upper=K> kk; // question for observation n
  array[N] int<lower=0, upper=1> y; // correctness for observation n
}
parameters {
  real delta; // mean student ability
  array[J] real alpha; // ability of student j - mean ability
  array[K] real beta; // difficulty of question k
}
model {
  alpha ~ normal(0, 1); // informative true prior
  beta ~ normal(0, 1); // informative true prior
  delta ~ normal(.75, 1); // informative true prior
  for (n in 1 : N) {
    y[n] ~ bernoulli_logit(alpha[jj[n]] - beta[kk[n]] + delta);
  }
}
