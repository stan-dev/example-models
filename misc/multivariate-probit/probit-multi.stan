functions {
  int sum(int[,] a) {
    int s;
    s <- 0;
    for (i in 1:size(a))
      for (j in 1:size(a[i]))
        s <- s + a[i,j];
    return s;
  }
}
data {
  int<lower=1> K;
  int<lower=1> D;
  int<lower=0> N;
  int<lower=0,upper=1> y[N,D];
  vector[K] x[N];
}
transformed data {
  int<lower=0> N_pos;
  int<lower=1,upper=N> n_pos[sum(y)];
  int<lower=1,upper=D> d_pos[size(n_pos)];
  int<lower=0> N_neg;
  int<lower=1,upper=N> n_neg[(N * D) - size(n_pos)];
  int<lower=1,upper=D> d_neg[size(n_neg)];

  N_pos <- size(n_pos);
  N_neg <- size(n_neg);
  {
    int i;
    int j;
    i <- 1;
    j <- 1;
    for (n in 1:N) {
      for (d in 1:D) {
        if (y[n,d] == 1) {
          n_pos[i] <- n;
          d_pos[i] <- d;
          i <- i + 1;
        } else {
          n_neg[j] <- n;
          d_neg[j] <- d;
          j <- j + 1;
        }
      }
    }
  }
}
parameters {
  matrix[D,K] beta;
  cholesky_factor_corr[D] L_Omega;
  vector<lower=0>[N_pos] z_pos;
  vector<upper=0>[N_neg] z_neg;
}
transformed parameters {
  vector[D] z[N];
  for (n in 1:N_pos)
    z[n_pos[n], d_pos[n]] <- z_pos[n];
  for (n in 1:N_neg)
    z[n_neg[n], d_neg[n]] <- z_neg[n];
}
model {
  L_Omega ~ lkj_corr_cholesky(4);
  to_vector(beta) ~ normal(0, 5);
  { 
    vector[D] beta_x[N];
    for (n in 1:N) 
      beta_x[n] <- beta * x[n];
    z ~ multi_normal_cholesky(beta_x, L_Omega);
  }
}
generated quantities {
  corr_matrix[D] Omega;
  Omega <- multiply_lower_tri_self_transpose(L_Omega);  
}
