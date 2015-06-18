## ---- irt-1pl-stan ----
data {
  int<lower=0> I;
  int<lower=0> J;
  int<lower=0,upper=1> y[I,J];
}
parameters {
  vector[I] b;
  vector[J] theta;
}
model {
  theta ~ normal(0, 1);
  b ~ normal(-1, 2);
  for (i in 1:I)
    y[i] ~ bernoulli_logit(theta - b[i]);
}
generated quantities {
  int<lower=0,upper=I> z_sim[100];
  real theta_sim[40];
  for (j in 1:40) {
    theta_sim[j] <- (j - 20) / 4.0;
    z_sim[j] <- 0;
    for (i in 1:I) {
      if (bernoulli_logit_rng(theta_sim[j] - b[i]))
        z_sim[j] <- z_sim[j] + 1;
    }
  }
}
