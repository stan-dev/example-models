functions {
  vector gpcm_probs(real theta, real mu, real alpha, vector beta) {
    vector[rows(beta) + 1] unsummed;
    vector[rows(beta) + 1] probs;
    unsummed <- append_row(rep_vector(0.0, 1), alpha*(theta + mu - beta));
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
  int<lower=0> y[N];             // response for n; y = 0, 1 ... m_i
  int<lower=1> K;                // # person covariates
  matrix[J,K] W;                 // person covariate matrix
}
transformed data {
  int r[N];                      // modified response; r = 1, 2, ... m_i + 1
  int m[I];                      // # parameters per item
  int pos[I];                    // first position in beta vector for item
  m <- rep_array(0, I);
  for(n in 1:N) {
    r[n] <- y[n] + 1;
    if(y[n] > m[ii[n]]) m[ii[n]] <- y[n];
  }
  pos[1] <- 1;
  for(i in 2:(I))
    pos[i] <- m[i-1] + pos[i-1];
}
parameters {
  vector<lower=0>[I] alpha;
  vector[sum(m)-1] beta_free;    // unconstrained item parameters
  vector[J] theta;
  vector[K] lambda;
}
transformed parameters {
  vector[sum(m)] beta;           // constrained item parameters
  beta <- append_row(beta_free, rep_vector(-1*sum(beta_free), 1));
}
model {
  vector[J] mu;
  mu <- W*lambda;
  theta ~ normal(0, 1);
  alpha ~ lognormal(1, 1);
  beta_free ~ normal(0, 5);
  for (n in 1:N)
    r[n] ~ categorical(gpcm_probs(theta[jj[n]], mu[jj[n]], alpha[ii[n]],
                                  segment(beta, pos[ii[n]], m[ii[n]])));
}
