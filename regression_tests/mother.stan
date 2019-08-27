// gen_test_mother1
// tests functions from current impelmentation:
//      generate_functions(prog.function_decl_defs_, o);
//      generate_member_var_decls_all(prog, o);
//      generate_constructor(prog, model_name, o);
//      generate_transform_inits_method(prog.parameter_decl_, o);

functions {

  int foo(int n);

  int foo(int n) {
    if (n == 0) return 1;
    return n * foo(n - 1);
  }


  real[] sho(real t,
             real[] y,
             real[] theta,
             data real[] x,
             data int[] x_int) ;

  real[] sho(real t,
             real[] y,
             real[] theta,
             data real[] x,
             data int[] x_int) {
    real dydt[2];
    dydt[1] = y[2];
    dydt[2] = -y[1] - theta[1] * y[2];
    return dydt;
  }

  real foo_bar0() {
    return 0.0;
  }

  real foo_bar1(real x) {
    return 1.0;
  }

  real foo_bar2(real x, real y) {
    return 2.0;
  }

  real foo_lpmf(int y, real lambda) {
    return 1.0;
  }

  real foo_lcdf(int y, real lambda) {
    return 1.0;
  }

  real foo_lccdf(int y, real lambda) {
    return 1.0;
  }

  void unit_normal_lp(real u) {
    increment_log_prob(normal_log(u,0,1));
    u ~ uniform(-100,100);
  }

  int foo_1(int a) {
    // direct while
    while (1) break;
    while (0) continue;

    // direct for
    for (i in 1:10) break;
    for (i in 1:10) continue;

    // in statement seq
    while (1) {
      int b;
      b = 5;
      break;
    }

    // if, else if, else body
    while (1) {
      if (0) break;
      else if (1) break;
      else break;
    }

    // nested while
    while (1) while (0) break;

    // nested for
    while (1) {
      for (i in 1:10) break;
    }

    // nested foreach (array)
    while (1) {
      int vs[2, 3];
      int z;
      for (v in vs) {
        z = 0;
        break;
      }
      for (v in vs) {
        z = 0;
        continue;
      }
      for (v in vs) {
        for (vv in v) {
          z = 0;
          break;
        }
        z = 1;
      }
    }

    // nested foreach (matrix)
    while (1) {
      real z;
      matrix[2,3] vs;
      for (v in vs) {
        z = 0;
        break;
      }
      for (v in vs) {
        z = 3.1;
        continue;
      }
    }

    // nested foreach (vector)
    while (1) {
      real z;
      vector[2] vs;
      for (v in vs) {
        z = 0;
        break;
      }
      for (v in vs) {
        z = 3.2;
        continue;
      }
    }

    // nested foreach (rowvector)
    while (1) {
      real z;
      row_vector[2] vs;
      for (v in vs) {
        z = 0;
        break;
      }
      for (v in vs) {
        z = 3.3;
        continue;
      }
    }

    // nested block
    while (1) {
      int b;
      b = 5;
      {
        int c;
        c = 6;
        break;
      }
    }

    return 0;
  }

  int foo_2(int a) {
    int vs[2];
    int y;
    for (v in vs) y = v;
    return 0;
  }

  real[] foo_3(real t, int n) {
    return rep_array(t,n);
  }

  real foo_lp(real x) {
    return x + get_lp();
  }

  void foo_4(real x) {
    reject("user-specified rejection", x);
  }

  real relative_diff(real x, real y, real max_, real min_) {
    real abs_diff;
    real avg_scale;
    abs_diff = fabs(x - y);
    avg_scale = (fabs(x) + fabs(y)) / 2;
    if ((abs_diff / avg_scale) > max_)
      reject("user-specified rejection, difference above ",max_," x:",x," y:",y);
    if ((abs_diff / avg_scale) < min_)
      reject("user-specified rejection, difference below ",min_," x:",x," y:",y);
    return abs_diff / avg_scale;
  }

  vector foo_5(vector shared_params, vector job_params,
             data real[] data_r, data int[] data_i) {
    return [1, 2, 3]';
  }

  real foo_five_args(real x1, real x2, real x3, real x4, real x5) {
    return x1;
  }
  real foo_five_args_lp(real x1, real x2, real x3, real x4, real x5, real x6) {
    return x1;
  }

  matrix covsqrt2corsqrt(matrix mat, int invert) {
    matrix[rows(mat), cols(mat)] o;
    o = mat;
    o[1] = o[2];
    o[3:4] = o[1:2];
    return o;
  }
  void f0(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    print("hi");
  }

  int f1(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a1;
  }

  int[] f2(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a2;
  }

  int[,] f3(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a3;
  }

  real f4(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a4;
  }

  real[] f5(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a5;
  }

  real[,] f6(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a6;
  }

  vector f7(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a7;
  }

  vector[] f8(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a8;
  }

  vector[,] f9(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a9;
  }

  matrix f10(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a10;
  }

  matrix[] f11(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a11;
  }

  matrix[,] f12(int a1, int[] a2, int[,] a3, real a4, real[] a5, real[,] a6, vector a7, vector[] a8, vector[,] a9, matrix a10, matrix[] a11, matrix[,] a12) {
    return a12;
  }
 void foo_6() {
    int a;
    real b;
    real c[20,30];
    matrix[40,50] ar_mat[60,70];
    ar_mat[1,1,1,1] = b;
  }
  matrix matfoo() {
    return [ [1,2,3,4,5,6,7,8,9,10]
             , [1,2,3,4,5,6,7,8,9,10]
             , [1,2,3,4,5,6,7,8,9,10] ];
  }
  vector vecfoo() {
    return [1,2,3,4,5,6,7,8,9,10]';
  }
  vector vecmufoo(real mu) {
    vector[10] l = mu * vecfoo();
    return l;
  }
  vector vecmubar(real mu) {
    vector[10] l = mu * [1,2,3,4,5,6,7,8,9,10]';
    return l[{1,2,3,4,5}];
  }
}
data {
  int<lower=0> N;
  int<lower=0> M;
  int<lower=0, upper=N*M> K;
  int<upper=N> d_int_1d_ar[N];
  int<upper=N> d_int_3d_ar[N,M,K];
  real<lower=-2.0, upper=2.0> J;
  real d_real_1d_ar[N];
  real d_real_3d_ar[N,M,K];
  vector[N] d_vec;
  vector[N] d_1d_vec[N];
  vector[N] d_3d_vec[N,M,K];
  row_vector[N] d_row_vec;
  row_vector[N] d_1d_row_vec[N];
  row_vector[N] d_3d_row_vec[N,M,K];
  matrix<lower=0,upper=1>[2,3] d_ar_mat[4,5];
  simplex[N] d_simplex;
  simplex[N] d_1d_simplex[N];
  simplex[N] d_3d_simplex[N,M,K];
  cholesky_factor_cov[5,4] d_cfcov_54;
  cholesky_factor_cov[3] d_cfcov_33;
  cholesky_factor_cov[3] d_cfcov_33_ar[K];
}
transformed data {
  int td_int;
  int td_1d[N];
  int td_1dk[M] = rep_array(1, M);
  int td_a = N;
  real td_b = N * J;
  real td_c = foo_bar1(td_b);
  matrix<lower=0,upper=1>[2,3] td_ar_mat[4,5];
  simplex[N] td_simplex;
  simplex[N] td_1d_simplex[N];
  simplex[N] td_3d_simplex[N,M,K];
  cholesky_factor_cov[5,5] td_cfcov_54; // TODO: Change to 5,4
  cholesky_factor_cov[3] td_cfcov_33;
  td_int = 1 || 2;
  td_int = 1 && 2;
  for (i in 1:2) {
    for (j in 1:3) {
      for (m in 1:4) {
        for (n in 1:5) {
          td_ar_mat[m, n, i, j] = 0.4;}}}}
  for (i in 1:N) {
    td_simplex[i] = 1.0 / N;
    for (n in 1:N) {
      td_1d_simplex[n, i] = 1.0 / N;
      for (m in 1:M) {
        for (k in 1:K) {
          td_3d_simplex[n, m, k, i] = 1.0 / N;}}}}
  for (i in 1:4) {
    for (j in 1:5) {
      matrix[2,3] l_mat = d_ar_mat[i,j];
      print("ar dim1: ",i, " ar dim2: ",j, " matrix: ", l_mat);
    }}
  td_cfcov_54 = diag_matrix(rep_vector(1, rows(td_cfcov_54)));
  td_cfcov_33 = diag_matrix(rep_vector(1, rows(td_cfcov_33)));
  {
    real z;
    row_vector[2] blocked_tdata_vs;
    for (v in blocked_tdata_vs) {
      z = 0;
    }
  }
  // some indexing tests for multi indices and slices
  td_1dk = td_1d[td_1dk];
  td_simplex = td_1d_simplex[1,:];
  td_simplex = td_1d_simplex[1,];
  td_simplex = td_1d_simplex[1,1:N];
}
parameters {
  real p_real;
  real<offset=1, multiplier=2> offset_multiplier[5];
  real<lower=0> p_real_1d_ar[N];
  real<lower=0> p_real_3d_ar[N,M,K];
  vector<lower=0>[N] p_vec;
  vector[N] p_1d_vec[N];
  vector[N] p_3d_vec[N,M,K];
  row_vector[N] p_row_vec;
  row_vector[N] p_1d_row_vec[N];
  row_vector[N] p_3d_row_vec[N,M,K];
  matrix<lower=0,upper=1>[2,3] p_ar_mat[4,5];
  simplex[N] p_simplex;
  simplex[N] p_1d_simplex[N];
  simplex[N] p_3d_simplex[N,M,K];
  cholesky_factor_cov[5,4] p_cfcov_54;
  cholesky_factor_cov[3] p_cfcov_33;
  cholesky_factor_cov[3] p_cfcov_33_ar[K];
}
transformed parameters {
  real<lower=0> tp_real_1d_ar[N];
  real<lower=0> tp_real_3d_ar[N,M,K];
  vector<upper=0>[N] tp_vec;
  vector[N] tp_1d_vec[N];
  vector[N] tp_3d_vec[N,M,K];
  row_vector[N] tp_row_vec;
  row_vector[N] tp_1d_row_vec[N];
  row_vector[N] tp_3d_row_vec[N,M,K];
  matrix<lower=0,upper=1>[2,3] tp_ar_mat[4,5];
  simplex[N] tp_simplex;
  simplex[N] tp_1d_simplex[N];
  simplex[N] tp_3d_simplex[N,M,K];
  cholesky_factor_cov[5,4] tp_cfcov_54;
  cholesky_factor_cov[3] tp_cfcov_33;
  cholesky_factor_cov[3] tp_cfcov_33_ar[K];

  tp_real_1d_ar = p_real_1d_ar;
  tp_real_3d_ar = p_real_3d_ar;
  tp_1d_vec = p_1d_vec;
  tp_3d_vec = p_3d_vec;

  tp_simplex = p_simplex;
  tp_1d_simplex = p_1d_simplex;
  tp_3d_simplex = p_3d_simplex;

  tp_cfcov_54 = p_cfcov_54;
  tp_cfcov_33 = p_cfcov_33;
  tp_cfcov_33_ar = p_cfcov_33_ar;

  for (i in 1:2) {
    for (j in 1:3) {
      for (m in 1:4) {
        for (n in 1:5) {
          tp_ar_mat[m, n, i, j] = 0.4;}}}}

  for (i in 1:N) tp_vec[i] = -1.0 * p_vec[i];
  tp_row_vec = tp_1d_vec[1]';
  tp_1d_row_vec = p_1d_row_vec;
  tp_3d_row_vec = p_3d_row_vec;
}
model {
  real r1 = foo_bar1(p_real);
  real r2 = foo_bar1(J);
  p_real ~ normal(0,1);
  offset_multiplier ~ normal(0, 1);

  to_vector(p_real_1d_ar) ~ normal(0, 1);
  for (n in 1:N) {
    to_vector(p_1d_vec[n]) ~ normal(0, 1);
    to_vector(p_1d_row_vec[n]) ~ normal(0, 1);
    to_vector(p_1d_simplex[n]) ~ normal(0, 1);
    for (m in 1:M) {
      for (k in 1:K) {
        to_vector(p_3d_vec[n, m, k]) ~ normal(d_3d_vec[n, m, k], 1);
        to_vector(p_3d_row_vec[n, m, k]) ~ normal(d_3d_row_vec[n, m, k], 1);
        to_vector(p_3d_simplex[n, m, k]) ~ normal(d_3d_simplex[n, m, k], 1);
        p_real_3d_ar[n, m, k] ~ normal(p_real_3d_ar[n, m, k], 1);
      }
    }
  }
  for (i in 1:4) {
    for (j in 1:5) {
      to_vector(p_ar_mat[i, j]) ~ normal(0, 1);
    }
  }
  for (k in 1:K) {
    to_vector(p_cfcov_33_ar[k]) ~ normal(0, 1);
  }
  to_vector(p_vec) ~ normal(d_vec, 1);
  to_vector(p_row_vec) ~ normal(0, 1);
  to_vector(p_simplex) ~ normal(0, 1);
  to_vector(p_cfcov_54) ~ normal(0, 1);
  to_vector(p_cfcov_33) ~ normal(0, 1);
}
generated quantities {
  real gq_r1 = foo_bar1(p_real);
  real gq_r2 = foo_bar1(J);
  real<lower = 0> gq_real_1d_ar[N];
  real<lower = 0> gq_real_3d_ar[N,M,K];
  vector<upper=1>[N] gq_vec;
  vector[N] gq_1d_vec[N];
  vector[N] gq_3d_vec[N,M,K];
  row_vector[N] gq_row_vec;
  row_vector[N] gq_1d_row_vec[N];
  row_vector[N] gq_3d_row_vec[N,M,K];
  matrix<lower=0,upper=1>[2,3] gq_ar_mat[4,5];
  simplex[N] gq_simplex;
  simplex[N] gq_1d_simplex[N];
  simplex[N] gq_3d_simplex[N,M,K];
  cholesky_factor_cov[5,4] gq_cfcov_54;
  cholesky_factor_cov[3] gq_cfcov_33;
  cholesky_factor_cov[3] gq_cfcov_33_ar[K];
  int indices[3] = {2, 3, 1};
  matrix[3, 4] indexing_mat[5];
  matrix[3, 4] idx_res1[3];
  matrix[3, 4] idx_res2[5];
  matrix[3, 3] idx_res3[3];
  matrix[3, 4] idx_res11[3];
  matrix[3, 4] idx_res21[5];
  matrix[3, 3] idx_res31[3];
  row_vector[4] idx_res4[3];
  vector[2] idx_res5[2];

  gq_real_1d_ar = p_1d_simplex[,1];
  gq_real_3d_ar = p_real_3d_ar;
  gq_1d_vec = p_1d_vec;
  gq_3d_vec = p_3d_vec;
  gq_row_vec = p_row_vec;
  gq_1d_row_vec = p_1d_row_vec;
  gq_3d_row_vec = p_3d_row_vec;

  gq_simplex = p_1d_simplex[1,1:N];
  gq_1d_simplex = p_1d_simplex;
  gq_3d_simplex = p_3d_simplex;

  gq_cfcov_54 = p_cfcov_54;
  gq_cfcov_33 = p_cfcov_33;
  gq_cfcov_33_ar = p_cfcov_33_ar;

  for (i in 1:2) {
    for (j in 1:3) {
      for (m in 1:4) {
        for (n in 1:5) {
          gq_ar_mat[m, n, i, j] = 0.4;}}}}

  for (i in 1:N) gq_vec[i] = -1.0 * p_vec[i];

  // A fun thing about Stan is that we can test syntactic sugar in Stan itself:
  for (i in 1:3)
    for (j in 1:4)
      for (k in 1:5)
        indexing_mat[k, i, j] = normal_rng(0, 1);

  // 2nd, 3rd, 1st indexing_matrix, 2nd, 3rd, 1st rows of each
  for (i in 1:size(indices))
    for (j in 1:size(indices))
      idx_res1[i, j] = indexing_mat[indices[i], indices[j]];

  idx_res11 = indexing_mat[indices, indices];
  if (indexing_mat[indices, indices][2,1,1] != idx_res1[2,1,1]) reject("indexing test 1 failed");

  //2nd, 3rd, 1st rows of every indexing_matrix
  for (i in 1:5)
    for (j in 1:size(indices))
      idx_res2[i, j] = indexing_mat[i, indices[j]];
  idx_res21 = indexing_mat[:, indices];
  //broken in stanc3
  if (indexing_mat[:, indices][2,1,1] != idx_res2[2,1,1]) reject("indexing test 2 failed");

  // (2nd, 3rd, 1st) indexing_matrices, all rows, 2nd, 3rd, 1st columns
  for (i in 1:size(indices))
    for (j in 1:3)
      for (k in 1:size(indices))
        idx_res3[i, j, k] = indexing_mat[indices[i], j, indices[k]];
  idx_res31 = indexing_mat[indices, :, indices];
  if (indexing_mat[indices, :, indices][2,1,1] != idx_res3[2,1,1]) reject("indexing test 3 failed");

  idx_res4 = indexing_mat[:3, 1, :];
  idx_res5 = indexing_mat[4:, 2:3, 1];
}
