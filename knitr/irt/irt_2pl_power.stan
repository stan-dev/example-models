## ---- irt-2pl-power-stan ----
data {
  int<lower=0> I;
  vector[I] a;
  vector[I] b;
}
model { 
}
generated quantities {
  int<lower=0,upper=I> z_sim[100];
  vector[100] theta_sim;
  for (j in 1:100) {
    theta_sim[j] <- (j - 50) / 10.0;
    z_sim[j] <- 0;
    for (i in 1:I)
      z_sim[j] <- z_sim[j] 
        + bernoulli_rng(inv_logit(a[i] * (theta_sim[j] - b[i])));
  }
}
