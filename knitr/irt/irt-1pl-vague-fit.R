## ---- fit-1pl-vague-stan ----
model_vague <- stan_model("irt_1pl_vague.stan");

mu_theta <- 0;  sigma_theta <- 100;
mu_b <- 0; sigma_b <- 100;
fit_vague <- sampling(model_vague, 
                 data=c("I", "J", "y", "mu_theta", "sigma_theta", "mu_b", "sigma_b"),
                 init=init_fun, refresh=2000, seed=1234);

## ---- fit-print-1pl-vague-stan ----
print(fit_vague, c(paste("b[", 1:10, "]"), paste("theta[", 1:10, "]"), "lp__"), 
      probs=c(0.10, 0.5, 0.90));

## ---- fit-print-1pl-vague-stan-one-fixed ----
mu_theta <- 0;  sigma_theta <- 1;
mu_b <- 0; sigma_b <- 100;
fit_vague2 <- sampling(model_vague, 
                  data=c("I", "J", "y", "mu_theta", "sigma_theta", "mu_b", "sigma_b"),
                  init=init_fun, refresh=2000, seed=1234);
print(fit_vague2, c(paste("b[", 1:10, "]"), paste("theta[", 1:10, "]"), "lp__"), 
      probs=c(0.10, 0.5, 0.90));
