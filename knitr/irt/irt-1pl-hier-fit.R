# ---- fit-1pl-hier-stan ----
library(rstan);
model_hier <- stan_model("irt_1pl_hier.stan");

fit_hier <-sampling(model_hier, data = c("I", "J", "y"),
                    refresh=2000, seed=1234)

options("width"=100);
print(fit_hier, c(paste("b[", 1:5, "]"), paste("theta[", 1:5, "]"),
                  "mu_b", "sigma_b", "sigma_theta", "lp__"), 
      probs=c(0.10, 0.5, 0.90));

print(b[1:10], digits=2);
print(theta[1:10], digits=2);
print("mu_b = 1;  sigma_b = 2;  sigma_theta=1");
