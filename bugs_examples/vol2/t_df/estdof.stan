//# estimated the degree of freedom (dof) from 
//# samples by sim.stan
// http://www.openbugs.net/Examples/t-df.html

//# estimate dof using continuous priors 

data {
  int<lower=0> N;
  array[N] real y;
}
parameters {
  // learning about the dof as a continuous quantity
  real<lower=2, upper=100> d;
}
model {
  y ~ student_t(d, 0, 1);
  // d ~ uniform(2, 100); 
}
