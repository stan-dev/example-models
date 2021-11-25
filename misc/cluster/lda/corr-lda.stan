data {
  int<lower=2> K;               // num topics
  int<lower=2> V;               // num words
  int<lower=1> M;               // num docs
  int<lower=0> corpus[M,V];     // word freq matrix, doc x word
  vector<lower=0>[V] beta;      // word prior
}
parameters {
  vector[K] mu;              // topic mean
  corr_matrix[K] Omega;      // correlation matrix
  vector<lower=0>[K] sigma;  // scales
  vector[K] eta[M];          // logit topic dist for doc m
  simplex[V] phi[K];         // word dist for topic k
}
transformed parameters {
  simplex[K] theta[M];       // simplex topic dist for doc m
  cov_matrix[K] Sigma;       // covariance matrix
  for (m in 1:M)
    theta[m] = softmax(eta[m]);
  for (m in 1:K) {
    Sigma[m,m] = sigma[m] * sigma[m] * Omega[m,m];
    for (n in (m+1):K) {
      Sigma[m,n] = sigma[m] * sigma[n] * Omega[m,n];
      Sigma[n,m] = Sigma[m,n];
    }
  } 
}
model {
  // priors
  for (k in 1:K)  
    phi[k] ~ dirichlet(beta);
  mu ~ normal(0,5);
  Omega ~ lkj_corr(2.0);
  sigma ~ cauchy(0,5);
  // topic distribution for docs
  for (m in 1:M)
    eta[m] ~ multi_normal(mu,Sigma);
  // token probabilities
  for (i in 1:M) {
    for (j in 1:V) {
      int count = corpus[i,j];
      real gamma[K];
      if (count > 0) {
        for (k in 1:K) {
          gamma[k] = (log(theta[i,k]) + log(phi[k,j]))*count;
        }
        increment_log_prob(log_sum_exp(gamma));  // likelihood
      }
    }
  }
}
