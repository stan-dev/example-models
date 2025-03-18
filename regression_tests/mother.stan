// gen_test_mother1
// tests functions from current impelmentation:
//      generate_functions(prog.function_decl_defs_, o);
//      generate_member_var_decls_all(prog, o);
//      generate_constructor(prog, model_name, o);
//      generate_transform_inits_method(prog.parameter_decl_, o);

functions {
  int foo(int n) {
    if (n == 0) {
      return 1;
    }
    return n * foo(n - 1);
  }
  
  array[] real sho(real t, array[] real y, array[] real theta,
                   data array[] real x, data array[] int x_int) {
    array[2] real dydt;
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
  
  real foo_rng(real mu, real sigma) {
    return normal_rng(mu, sigma);
  }
  
  void unit_normal_lp(real u) {
    target += normal_lpdf(u | 0, 1);
    u ~ uniform(-100, 100);
  }
  
  int foo_1(int a) {
    // direct while
    while (1) {
      break;
    }
    while (0) {
      continue;
    }
    
    // direct for
    for (i in 1 : 10) {
      break;
    }
    for (i in 1 : 10) {
      continue;
    }
    
    // in statement seq
    while (1) {
      int b;
      b = 5;
      break;
    }
    
    // if, else if, else body
    while (1) {
      if (0) {
        break;
      } else if (1) {
        break;
      } else {
        break;
      }
    }
    
    // nested while
    while (1) {
      while (0) {
        break;
      }
    }
    
    // nested for
    while (1) {
      for (i in 1 : 10) {
        break;
      }
    }
    
    // nested foreach (array)
    while (1) {
      array[2, 3] int vs;
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
      matrix[2, 3] vs;
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
    array[2] int vs;
    int y;
    for (v in vs) 
      y = v;
    return 0;
  }
  
  array[] real foo_3(real t, int n) {
    return rep_array(t, n);
  }
  
  real foo_lp(real x) {
    return x + target();
  }
  
  void foo_4(real x) {
    reject("user-specified rejection", x);
  }
  
  real relative_diff(real x, real y, real max_, real min_) {
    real abs_diff;
    real avg_scale;
    abs_diff = abs(x - y);
    avg_scale = (abs(x) + abs(y)) / 2;
    if ((abs_diff / avg_scale) > max_) {
      reject("user-specified rejection, difference above ", max_, " x:", x,
             " y:", y);
    }
    if ((abs_diff / avg_scale) < min_) {
      reject("user-specified rejection, difference below ", min_, " x:", x,
             " y:", y);
    }
    return abs_diff / avg_scale;
  }
  
  vector foo_5(vector shared_params, vector job_params,
               data array[] real data_r, data array[] int data_i) {
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
    o[3 : 4] = o[1 : 2];
    return o;
  }
  void f0(int a1, array[] int a2, array[,] int a3, real a4, array[] real a5,
          array[,] real a6, vector a7, array[] vector a8, array[,] vector a9,
          matrix a10, array[] matrix a11, array[,] matrix a12) {
    print("hi");
  }
  
  int f1(int a1, array[] int a2, array[,] int a3, real a4, array[] real a5,
         array[,] real a6, vector a7, array[] vector a8, array[,] vector a9,
         matrix a10, array[] matrix a11, array[,] matrix a12) {
    return a1;
  }
  
  array[] int f2(int a1, array[] int a2, array[,] int a3, real a4,
                 array[] real a5, array[,] real a6, vector a7,
                 array[] vector a8, array[,] vector a9, matrix a10,
                 array[] matrix a11, array[,] matrix a12) {
    return a2;
  }
  
  array[,] int f3(int a1, array[] int a2, array[,] int a3, real a4,
                  array[] real a5, array[,] real a6, vector a7,
                  array[] vector a8, array[,] vector a9, matrix a10,
                  array[] matrix a11, array[,] matrix a12) {
    return a3;
  }
  
  real f4(int a1, array[] int a2, array[,] int a3, real a4, array[] real a5,
          array[,] real a6, vector a7, array[] vector a8, array[,] vector a9,
          matrix a10, array[] matrix a11, array[,] matrix a12) {
    return a4;
  }
  
  array[] real f5(int a1, array[] int a2, array[,] int a3, real a4,
                  array[] real a5, array[,] real a6, vector a7,
                  array[] vector a8, array[,] vector a9, matrix a10,
                  array[] matrix a11, array[,] matrix a12) {
    return a5;
  }
  
  array[,] real f6(int a1, array[] int a2, array[,] int a3, real a4,
                   array[] real a5, array[,] real a6, vector a7,
                   array[] vector a8, array[,] vector a9, matrix a10,
                   array[] matrix a11, array[,] matrix a12) {
    return a6;
  }
  
  vector f7(int a1, array[] int a2, array[,] int a3, real a4,
            array[] real a5, array[,] real a6, vector a7, array[] vector a8,
            array[,] vector a9, matrix a10, array[] matrix a11,
            array[,] matrix a12) {
    return a7;
  }
  
  array[] vector f8(int a1, array[] int a2, array[,] int a3, real a4,
                    array[] real a5, array[,] real a6, vector a7,
                    array[] vector a8, array[,] vector a9, matrix a10,
                    array[] matrix a11, array[,] matrix a12) {
    return a8;
  }
  
  array[,] vector f9(int a1, array[] int a2, array[,] int a3, real a4,
                     array[] real a5, array[,] real a6, vector a7,
                     array[] vector a8, array[,] vector a9, matrix a10,
                     array[] matrix a11, array[,] matrix a12) {
    return a9;
  }
  
  matrix f10(int a1, array[] int a2, array[,] int a3, real a4,
             array[] real a5, array[,] real a6, vector a7, array[] vector a8,
             array[,] vector a9, matrix a10, array[] matrix a11,
             array[,] matrix a12) {
    return a10;
  }
  
  array[] matrix f11(int a1, array[] int a2, array[,] int a3, real a4,
                     array[] real a5, array[,] real a6, vector a7,
                     array[] vector a8, array[,] vector a9, matrix a10,
                     array[] matrix a11, array[,] matrix a12) {
    return a11;
  }
  
  array[,] matrix f12(int a1, array[] int a2, array[,] int a3, real a4,
                      array[] real a5, array[,] real a6, vector a7,
                      array[] vector a8, array[,] vector a9, matrix a10,
                      array[] matrix a11, array[,] matrix a12) {
    return a12;
  }
  void foo_6() {
    int a;
    real b;
    array[20, 30] real c;
    array[60, 70] matrix[40, 50] ar_mat;
    ar_mat[1, 1, 1, 1] = b;
  }
  matrix matfoo() {
    return [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]];
  }
  vector vecfoo() {
    return [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]';
  }
  vector vecmufoo(real mu) {
    vector[10] l = mu * vecfoo();
    return l;
  }
  vector vecmubar(real mu) {
    vector[10] l = mu * [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]';
    return l[{1, 2, 3, 4, 5}];
  }
  vector algebra_system(vector x, vector y, array[] real dat,
                        array[] int dat_int) {
    vector[2] f_x;
    f_x[1] = x[1] - y[1];
    f_x[2] = x[2] - y[2];
    return f_x;
  }
  
  vector binomialf(vector phi, vector theta, data array[] real x_r,
                   data array[] int x_i) {
    vector[1] lpmf;
    lpmf[1] = 0.0;
    return lpmf;
  }
}
data {
  int<lower=0> N;
  int<lower=0> M;
  int<lower=0, upper=N * M> K;
  array[N] int<upper=N> d_int_1d_ar;
  array[N, M, K] int<upper=N> d_int_3d_ar;
  real<lower=-2.0, upper=2.0> J;
  array[N] real d_real_1d_ar;
  array[N, M, K] real d_real_3d_ar;
  vector[N] d_vec;
  array[N] vector[N] d_1d_vec;
  array[N, M, K] vector[N] d_3d_vec;
  row_vector[N] d_row_vec;
  array[N] row_vector[N] d_1d_row_vec;
  array[N, M, K] row_vector[N] d_3d_row_vec;
  array[4, 5] matrix<lower=0, upper=1>[2, 3] d_ar_mat;
  simplex[N] d_simplex;
  array[N] simplex[N] d_1d_simplex;
  array[N, M, K] simplex[N] d_3d_simplex;
  cholesky_factor_cov[5, 4] d_cfcov_54;
  cholesky_factor_cov[3] d_cfcov_33;
  array[K] cholesky_factor_cov[3] d_cfcov_33_ar;
}
transformed data {
  int td_int;
  array[N] int td_1d;
  array[M] int td_1dk = rep_array(1, M);
  int td_a = N;
  real td_b = N * J;
  real td_c = foo_bar1(td_b);
  real td_d = foo_rng(td_b, td_c);
  array[4, 5] matrix<lower=0, upper=1>[2, 3] td_ar_mat;
  simplex[N] td_simplex;
  array[N] simplex[N] td_1d_simplex;
  array[N, M, K] simplex[N] td_3d_simplex;
  cholesky_factor_cov[5, 5] td_cfcov_54; // TODO: Change to 5,4
  cholesky_factor_cov[3] td_cfcov_33;
  vector[2] x = [0, 0]';
  vector[2] y = [1, 1]';
  array[0] real dat;
  array[0] int dat_int;
  array[0, 0] real x_r;
  array[0, 0] int x_i;
  td_int = 1 || 2;
  td_int = 1 && 2;
  for (i in 1 : 2) {
    for (j in 1 : 3) {
      for (m in 1 : 4) {
        for (n in 1 : 5) {
          td_ar_mat[m, n, i, j] = 0.4;
        }
      }
    }
  }
  for (i in 1 : N) {
    td_simplex[i] = 1.0 / N;
    for (n in 1 : N) {
      td_1d_simplex[n, i] = 1.0 / N;
      for (m in 1 : M) {
        for (k in 1 : K) {
          td_3d_simplex[n, m, k, i] = 1.0 / N;
        }
      }
    }
  }
  for (i in 1 : 4) {
    for (j in 1 : 5) {
      matrix[2, 3] l_mat = d_ar_mat[i, j];
      print("ar dim1: ", i, " ar dim2: ", j, " matrix: ", l_mat);
    }
  }
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
  td_simplex = td_1d_simplex[1,  : ];
  td_simplex = td_1d_simplex[1,  : ];
  td_simplex = td_1d_simplex[1, 1 : N];
  
  array[2, 2] int arr_mul_ind;
  arr_mul_ind[1, 1 : 2] = {1, 1};
  
  array[2] real x_mul_ind = {1, 2};
}
parameters {
  real p_real;
  array[5] real<offset=1, multiplier=2> offset_multiplier;
  array[N] real<lower=0> p_real_1d_ar;
  array[N, M, K] real<lower=0> p_real_3d_ar;
  vector<lower=0>[N] p_vec;
  array[N] vector[N] p_1d_vec;
  array[N, M, K] vector[N] p_3d_vec;
  row_vector[N] p_row_vec;
  array[N] row_vector[N] p_1d_row_vec;
  array[N, M, K] row_vector[N] p_3d_row_vec;
  array[4, 5] matrix<lower=0, upper=1>[2, 3] p_ar_mat;
  simplex[N] p_simplex;
  array[N] simplex[N] p_1d_simplex;
  array[N, M, K] simplex[N] p_3d_simplex;
  cholesky_factor_cov[5, 4] p_cfcov_54;
  cholesky_factor_cov[3] p_cfcov_33;
  array[K] cholesky_factor_cov[3] p_cfcov_33_ar;
  vector[2] x_p;
  vector[2] y_p;
}
transformed parameters {
  array[N] real<lower=0> tp_real_1d_ar;
  array[N, M, K] real<lower=0> tp_real_3d_ar;
  vector<upper=0>[N] tp_vec;
  array[N] vector[N] tp_1d_vec;
  array[N, M, K] vector[N] tp_3d_vec;
  row_vector[N] tp_row_vec;
  array[N] row_vector[N] tp_1d_row_vec;
  array[N, M, K] row_vector[N] tp_3d_row_vec;
  array[4, 5] matrix<lower=0, upper=1>[2, 3] tp_ar_mat;
  simplex[N] tp_simplex;
  array[N] simplex[N] tp_1d_simplex;
  array[N, M, K] simplex[N] tp_3d_simplex;
  cholesky_factor_cov[5, 4] tp_cfcov_54;
  cholesky_factor_cov[3] tp_cfcov_33;
  array[K] cholesky_factor_cov[3] tp_cfcov_33_ar;
  vector[2] theta_p;
  
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
  
  for (i in 1 : 2) {
    for (j in 1 : 3) {
      for (m in 1 : 4) {
        for (n in 1 : 5) {
          tp_ar_mat[m, n, i, j] = 0.4;
        }
      }
    }
  }
  
  for (i in 1 : N) {
    tp_vec[i] = -1.0 * p_vec[i];
  }
  tp_row_vec = tp_1d_vec[1]';
  tp_1d_row_vec = p_1d_row_vec;
  tp_3d_row_vec = p_3d_row_vec;
  
  theta_p = algebra_solver(algebra_system, x, y, dat, dat_int);
  theta_p = algebra_solver(algebra_system, x, y, dat, dat_int, 0.01, 0.01,
                           10);
  theta_p = algebra_solver(algebra_system, x, y_p, dat, dat_int, 0.01, 0.01,
                           10);
  theta_p = algebra_solver(algebra_system, x_p, y, dat, dat_int);
  theta_p = algebra_solver(algebra_system, x_p, y, dat, dat_int, 0.01, 0.01,
                           10);
  theta_p = algebra_solver(algebra_system, x_p, y_p, dat, dat_int);
  theta_p = algebra_solver(algebra_system, x_p, y_p, dat, dat_int, 0.01,
                           0.01, 10);
}
model {
  vector[0] tmp;
  array[0] vector[0] tmp2;
  real r1 = foo_bar1(p_real);
  real r2 = foo_bar1(J);
  unit_normal_lp(p_real);
  p_real ~ normal(0, 1);
  offset_multiplier ~ normal(0, 1);
  
  to_vector(p_real_1d_ar) ~ normal(0, 1);
  for (n in 1 : N) {
    to_vector(p_1d_vec[n]) ~ normal(0, 1);
    to_vector(p_1d_row_vec[n]) ~ normal(0, 1);
    to_vector(p_1d_simplex[n]) ~ normal(0, 1);
    to_vector(to_matrix(p_real_3d_ar[n])) ~ normal(0, 1);

    for (m in 1 : M) {
      for (k in 1 : K) {
        to_vector(p_3d_vec[n, m, k]) ~ normal(d_3d_vec[n, m, k], 1);
        to_vector(p_3d_row_vec[n, m, k]) ~ normal(d_3d_row_vec[n, m, k], 1);
        to_vector(p_3d_simplex[n, m, k]) ~ normal(d_3d_simplex[n, m, k], 1);
        p_real_3d_ar[n, m, k] ~ normal(p_real_3d_ar[n, m, k], 1);
      }
    }
  }
  for (i in 1 : 4) {
    for (j in 1 : 5) {
      to_vector(p_ar_mat[i, j]) ~ normal(0, 1);
    }
  }
  for (k in 1 : K) {
    to_vector(p_cfcov_33_ar[k]) ~ normal(0, 1);
  }
  to_vector(p_vec) ~ normal(d_vec, 1);
  to_vector(p_row_vec) ~ normal(0, 1);
  to_vector(p_simplex) ~ normal(0, 1);
  to_vector(p_cfcov_54) ~ normal(0, 1);
  to_vector(p_cfcov_33) ~ normal(0, 1);
  
  target += map_rect(binomialf, tmp, tmp2, x_r, x_i);
}
generated quantities {
  real gq_r1 = foo_bar1(p_real);
  real gq_r2 = foo_bar1(J);
  array[N] real<lower=0> gq_real_1d_ar;
  array[N, M, K] real<lower=0> gq_real_3d_ar;
  vector<upper=1>[N] gq_vec;
  array[N] vector[N] gq_1d_vec;
  array[N, M, K] vector[N] gq_3d_vec;
  row_vector[N] gq_row_vec;
  array[N] row_vector[N] gq_1d_row_vec;
  array[N, M, K] row_vector[N] gq_3d_row_vec;
  array[4, 5] matrix<lower=0, upper=1>[2, 3] gq_ar_mat;
  simplex[N] gq_simplex;
  array[N] simplex[N] gq_1d_simplex;
  array[N, M, K] simplex[N] gq_3d_simplex;
  cholesky_factor_cov[5, 4] gq_cfcov_54;
  cholesky_factor_cov[3] gq_cfcov_33;
  array[K] cholesky_factor_cov[3] gq_cfcov_33_ar;
  array[3] int indices = {2, 3, 1};
  array[5] matrix[3, 4] indexing_mat;
  array[3] matrix[3, 4] idx_res1;
  array[5] matrix[3, 4] idx_res2;
  array[3] matrix[3, 3] idx_res3;
  array[3] matrix[3, 4] idx_res11;
  array[5] matrix[3, 4] idx_res21;
  array[3] matrix[3, 3] idx_res31;
  array[3] row_vector[4] idx_res4;
  array[2] vector[2] idx_res5;
  
  gq_real_1d_ar = p_1d_simplex[ : , 1];
  gq_real_3d_ar = p_real_3d_ar;
  gq_1d_vec = p_1d_vec;
  gq_3d_vec = p_3d_vec;
  gq_row_vec = p_row_vec;
  gq_1d_row_vec = p_1d_row_vec;
  gq_3d_row_vec = p_3d_row_vec;
  
  gq_simplex = p_1d_simplex[1, 1 : N];
  gq_1d_simplex = p_1d_simplex;
  gq_3d_simplex = p_3d_simplex;
  
  gq_cfcov_54 = p_cfcov_54;
  gq_cfcov_33 = p_cfcov_33;
  gq_cfcov_33_ar = p_cfcov_33_ar;
  
  for (i in 1 : 2) {
    for (j in 1 : 3) {
      for (m in 1 : 4) {
        for (n in 1 : 5) {
          gq_ar_mat[m, n, i, j] = 0.4;
        }
      }
    }
  }
  
  for (i in 1 : N) {
    gq_vec[i] = -1.0 * p_vec[i];
  }
  
  // A fun thing about Stan is that we can test syntactic sugar in Stan itself:
  for (i in 1 : 3) {
    for (j in 1 : 4) {
      for (k in 1 : 5) {
        indexing_mat[k, i, j] = normal_rng(0, 1);
      }
    }
  }
  
  // 2nd, 3rd, 1st indexing_matrix, 2nd, 3rd, 1st rows of each
  for (i in 1 : size(indices)) {
    for (j in 1 : size(indices)) {
      idx_res1[i, j] = indexing_mat[indices[i], indices[j]];
    }
  }
  
  idx_res11 = indexing_mat[indices, indices];
  if (indexing_mat[indices, indices][2, 1, 1] != idx_res1[2, 1, 1]) {
    reject("indexing test 1 failed");
  }
  
  //2nd, 3rd, 1st rows of every indexing_matrix
  for (i in 1 : 5) {
    for (j in 1 : size(indices)) {
      idx_res2[i, j] = indexing_mat[i, indices[j]];
    }
  }
  idx_res21 = indexing_mat[ : , indices];
  //broken in stanc3
  if (indexing_mat[ : , indices][2, 1, 1] != idx_res2[2, 1, 1]) {
    reject("indexing test 2 failed");
  }
  
  // (2nd, 3rd, 1st) indexing_matrices, all rows, 2nd, 3rd, 1st columns
  for (i in 1 : size(indices)) {
    for (j in 1 : 3) {
      for (k in 1 : size(indices)) {
        idx_res3[i, j, k] = indexing_mat[indices[i], j, indices[k]];
      }
    }
  }
  idx_res31 = indexing_mat[indices,  : , indices];
  if (indexing_mat[indices,  : , indices][2, 1, 1] != idx_res3[2, 1, 1]) {
    reject("indexing test 3 failed");
  }
  
  idx_res4 = indexing_mat[ : 3, 1,  : ];
  idx_res5 = indexing_mat[4 : , 2 : 3, 1];
}
