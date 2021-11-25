data {
  int<lower=2> K;               // num topics
  int<lower=2> V;               // num words
  int<lower=1> M;               // num docs
  int<lower=0> corpus[M,V];     // word freq matrix, doc x word
  vector<lower=0>[K] alpha;     // topic prior
  vector<lower=0>[V] beta;      // word prior
}
parameters {
  simplex[K] theta[M];   // topic dist for doc m
  simplex[V] phi[K];     // word dist for topic k
}
model {
  for (m in 1:M)  
    theta[m] ~ dirichlet(alpha);  // prior
  for (k in 1:K)  
    phi[k] ~ dirichlet(beta);     // prior
  for (i in 1:M) {
    for (j in 1:V) {
      int count = corpus[i,j];
      real gamma[K];
      if (count > 0) {
        for (k in 1:K) {
          gamma[k] <- (log(theta[i,k]) + log(phi[k,j]))*count;
        }
        increment_log_prob(log_sum_exp(gamma));  // likelihood
      }
    }
  }
}
