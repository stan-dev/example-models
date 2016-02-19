data {
  int<lower=1> I;               // # questions
  int<lower=1> J;               // # persons
  int<lower=1> N;               // # observations
  int<lower=1, upper=I> ii[N];  // question for n
  int<lower=1, upper=J> jj[N];  // person for n
  int<lower=0, upper=1> y[N];   // correctness for n
  int<lower=1> K;               // # person covariates
  matrix[J,K] W;                // person covariate matrix
}
parameters {
  vector<lower=0>[I] alpha;
  vector[I-1] beta_free;
  vector[J] theta;
  vector[K] lambda;
}
transformed parameters {
  vector[I] beta;
  beta <- append_row(beta_free, rep_vector(-1*sum(beta_free), 1));
}
model {
  vector[J] mu;
  mu <- W*lambda;
  alpha ~ lognormal(1, 1);
  beta_free ~ normal(0, 5);
  theta ~ normal(0, 1);
  y ~ bernoulli_logit(alpha[ii].*theta[jj] + mu[jj] - beta[ii]);
}
