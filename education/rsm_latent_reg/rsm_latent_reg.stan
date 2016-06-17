functions {
  vector rsm_probs(real theta, real beta, vector kappa) {
    vector[rows(kappa) + 1] unsummed;
    vector[rows(kappa) + 1] probs;
    unsummed <- append_row(rep_vector(0, 1), theta - beta - kappa);
    probs <- softmax(cumulative_sum(unsummed));
    return probs;
  }
}
data {
  int<lower=1> I;                // # items
  int<lower=1> J;                // # persons
  int<lower=1> N;                // # responses
  int<lower=1,upper=I> ii[N];    // i for n
  int<lower=1,upper=J> jj[N];    // j for n
  int<lower=0> y[N];             // response for n; y in {0 ... m_i}
  int<lower=1> K;                // # person covariates
  matrix[J,K] W;                 // person covariate matrix
}
transformed data {
  int r[N];                      // modified response; r in {1 ... m_i + 1}
  int m;                         // # steps
  m <- max(y);
  for(n in 1:N)
    r[n] <- y[n] + 1;
}
parameters {
  vector[I-1] beta_free;         // unconstrained item parameters
  vector[m-1] kappa_free;        // unconstrained step parameters
  vector[J] theta;
  real<lower=0> sigma;
  vector[K] lambda;
}
transformed parameters {
  vector[I] beta;                // all item parameters
  vector[m] kappa;               // all step parameters
  beta <- append_row(beta_free, rep_vector(-1*sum(beta_free), 1));
  kappa <- append_row(kappa_free, rep_vector(-1*sum(kappa_free), 1));
}
model {
  theta ~ normal(W*lambda, sigma);
  sigma ~ exponential(.1);
  beta_free ~ normal(0, 5);
  kappa_free ~ normal(0, 5);
  for (n in 1:N)
    r[n] ~ categorical(rsm_probs(theta[jj[n]], beta[ii[n]], kappa));
}
