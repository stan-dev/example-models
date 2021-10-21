functions {
  int sum_2darray(array[,] int a) {
    int s = 0;
    for (i in 1 : size(a)) {
      s = s + sum(a[i]);
    }
    return s;
  }
}
data {
  int<lower=1> K;
  int<lower=1> D;
  int<lower=0> N;
  array[N, D] int<lower=0, upper=1> y;
  array[N] vector[K] x;
}
transformed data {
  int<lower=0> N_pos;
  array[sum_2darray(y)] int<lower=1, upper=N> n_pos;
  array[size(n_pos)] int<lower=1, upper=D> d_pos;
  int<lower=0> N_neg;
  array[(N * D) - size(n_pos)] int<lower=1, upper=N> n_neg;
  array[size(n_neg)] int<lower=1, upper=D> d_neg;
  
  N_pos = size(n_pos);
  N_neg = size(n_neg);
  {
    int i = 1;
    int j = 1;
    for (n in 1 : N) {
      for (d in 1 : D) {
        if (y[n, d] == 1) {
          n_pos[i] = n;
          d_pos[i] = d;
          i = i + 1;
        } else {
          n_neg[j] = n;
          d_neg[j] = d;
          j = j + 1;
        }
      }
    }
  }
}
parameters {
  matrix[D, K] beta;
  cholesky_factor_corr[D] L_Omega;
  vector<lower=0>[N_pos] z_pos;
  vector<upper=0>[N_neg] z_neg;
}
transformed parameters {
  array[N] vector[D] z;
  for (n in 1 : N_pos) {
    z[n_pos[n], d_pos[n]] = z_pos[n];
  }
  for (n in 1 : N_neg) {
    z[n_neg[n], d_neg[n]] = z_neg[n];
  }
}
model {
  L_Omega ~ lkj_corr_cholesky(4);
  to_vector(beta) ~ normal(0, 5);
  {
    array[N] vector[D] beta_x;
    for (n in 1 : N) {
      beta_x[n] = beta * x[n];
    }
    z ~ multi_normal_cholesky(beta_x, L_Omega);
  }
}
generated quantities {
  corr_matrix[D] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
