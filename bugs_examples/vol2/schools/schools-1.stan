/**
 * Schools: ranking school examination resutls using 
 * http://www.openbugs.info/Examples/Schools.html
 */
data {
  int<lower=0> N; 
  int<lower=0> M; 
  vector[N] LRT;
  int school[N]; 
  int School_denom[N, 3]; 
  int School_gender[N, 2]; 
  int VR[N, 2]; 
  real Y[N]; 
  int Gender[N]; 
  cov_matrix[3] R; 
} 
parameters {
  real beta[8]; 
  real<lower=0> theta; 
  real phi; 
  vector[3] alpha[M];
  vector[3] gamma;
  cov_matrix[3] Omega;
} 
model {
  real Ymu[N]; 
  for(p in 1:N)
    Ymu[p] <- alpha[school[p],1] 
      + alpha[school[p],2] * LRT[p] 
      + alpha[school[p],3] * VR[p,1] 
      + beta[1] * LRT[p] * LRT[p] 
      + beta[2] * VR[p,2] 
      + beta[3] * Gender[p] 
      + beta[4] * School_gender[p,1]
      + beta[5] * School_gender[p,2]
      + beta[6] * School_denom[p,1] 
      + beta[7] * School_denom[p,2] 
      + beta[8] * School_denom[p,3];

  Y ~ normal(Ymu,  exp(-0.5 * (theta + phi * LRT)));

  # Priors for fixed effects:
  beta ~ normal(0, 5);
  theta ~ normal(0, 5);
  phi ~ normal(0, 5);

  # Priors for random coefficients:
  alpha ~ multi_normal_prec(gamma, Omega);

  # Hyper-priors:
  gamma ~ normal(0,5);
  Omega ~ wishart(3, R); 
}

