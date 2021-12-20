/*
 * Endo: conditional inference in case-contrl studies 
 * http://www.openbugs.net/Examples/Endo.html
 *
 * In this example, three methods of different 
 * model specifications are used for one equivalent
 * model. This is method 2. 

 * FIXME: using the multinomial specification  
 */
data {
  int n10;
  int n01;
  int n11;
  int I;
  int J;
}
transformed data {
  array[I, 2] int<lower=0> Y;
  array[I, 2] int<lower=0> est;
  for (i in 1 : I) {
    Y[i, 1] = 1;
    Y[i, 2] = 0;
  }
  for (i in 1 : n10) {
    est[i, 1] = 1;
    est[i, 2] = 0;
  }
  for (i in (n10 + 1) : (n10 + n01)) {
    est[i, 1] = 0;
    est[i, 2] = 1;
  }
  for (i in (n10 + n01 + 1) : (n10 + n01 + n11)) {
    est[i, 1] = 1;
    est[i, 2] = 1;
  }
  for (i in (n10 + n01 + n11 + 1) : I) {
    est[i, 1] = 0;
    est[i, 2] = 0;
  }
}
parameters {
  real beta;
}
model {
  array[I, 2] real p;
  
  // METHOD 2 - conditional likelihoods
  beta ~ normal(0, 1000);
  
  for (i in 1 : I) {
    p[i, 1] = exp(beta * est[i, 1]);
    p[i, 2] = exp(beta * est[i, 2]);
    p[i, 1] = p[i, 1] / (p[i, 1] + p[i, 2]);
    p[i, 2] = 1 - p[i, 1];
    
    // using the multinomial log-pmf explicitly 
    target += log(p[i, 1]) * Y[i, 1];
    target += log(p[i, 2]) * Y[i, 2];
  }
}
