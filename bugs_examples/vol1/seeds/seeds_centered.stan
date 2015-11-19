#  http://www.openbugs.net/Examples/Seeds.html
# tau replaced by sigma direct estimation
# using narrower semi-informative priors
# Centered coefficients as in BUGS manual chapter 3.1
data {
    int<lower=0> I;
    int<lower=0> n[I];
    int<lower=0> N[I];
    vector[I] x1; # seed type
    vector[I] x2; # root extract
}

transformed data {
    vector[I] x1x2;
    x1x2 <- x1 .* x2;
}
parameters {
    real alpha0;
    real alpha1;
    real alpha12;
    real alpha2;
    vector[I] c;
    real<lower=0> sigma;
}

transformed parameters{
  vector[I] b;
  b <- c - mean(c);
}

model {
   alpha0 ~  normal(0.0, 1.0); # Narrower priors
   alpha1 ~  normal(0.0, 1.0);
   alpha2 ~  normal(0.0, 1.0);
   alpha12 ~ normal(0.0, 1.0);
   sigma ~   cauchy(0,1);

   c ~ normal(0.0, sigma);
   n ~ binomial_logit(N, alpha0 + alpha1 * x1 + alpha2 * x2 + alpha12 * x1x2 + b);
}
