data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  vector[N] z;
}
transformed data {
  vector[N] y_std;
  vector[N] x_std;
  vector[N] z_std;

  y_std <- (y - mean(y)) / sd(y);
  x_std <- (x - mean(x)) / sd(x);
  z_std <- (z - mean(z)) / sd(z);
}
parameters {
  real alpha;
  real beta;
  real gamma;
  real delta;
  vector<lower=0>[2] sigma;
}
transformed parameters {
  vector[2] log_joint;

  log_joint[1] 
    <- normal_log(alpha,0,10)
    + normal_log(beta,0,5)
    + cauchy_log(sigma[1],0,5)
    + normal_log(y_std, alpha + beta * x_std, sigma[1]);
  log_joint[2]
    <- normal_log(gamma,0,10)
    + normal_log(delta,0,5)
    + cauchy_log(sigma[2],0,5)
    + normal_log(y_std, gamma + delta * z_std, sigma[2]);
}
model {
  increment_log_prob(log_joint[1]);
  increment_log_prob(log_joint[2]);
}
generated quantities {
  real lambda;
  lambda <- softmax(log_joint)[1];
}
