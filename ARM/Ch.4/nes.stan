data {
  int<lower=0> N;
  vector[N] partyid7;
  vector[N] real_ideo;
  vector[N] race_adj;
  vector[N] educ1;
  vector[N] gender;
  vector[N] income;
  int age_discrete[N];
}
transformed data {
  vector[N] age30_44;        // age as factor
  vector[N] age45_64;
  vector[N] age65up;
  matrix[N,8] x;

  for (n in 1:N) {
    age30_44[n] = age_discrete[n] == 2;
    age45_64[n] = age_discrete[n] == 3;
    age65up[n]  = age_discrete[n] == 4;
  }
  
  x = [real_ideo', race_adj', age30_44', age45_64',
       age65up', educ1', gender', income']';
}
parameters {
  real alpha;
  vector[8] beta;
  real<lower=0> sigma;
}
model {
  partyid7 ~ normal_id_glm(x, alpha, beta, sigma);
}
