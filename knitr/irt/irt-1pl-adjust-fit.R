## ---- fit-1pl-stan-adjust ----
model_adj <- stan_model("irt_1pl_adjust.stan");
fit_adj <- sampling(model_adj, data=c("I", "J", "y"), refresh=2000, seed=1234);

## ---- fit-print-1pl-stan-adjust ----
print(fit_adj, c(paste("theta_raw[", 1:5, "]"),
                 paste("theta[", 1:5, "]"),
                 paste("b_raw[", 1:5, "]"),
                 paste("b[", 1:5, "]"),
                 "lp__"), 
      probs=c(0.10, 0.5, 0.90));
