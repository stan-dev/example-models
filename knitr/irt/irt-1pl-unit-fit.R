## ---- fit-1pl-stan-unit ----
model_unit <- stan_model("irt_1pl_unit.stan");
fit_unit <- sampling(model_adj, data=c("I", "J", "y"), refresh=2000, seed=1234);

## ---- fit-print-1pl-stan-unit ----
print(fit_unit, c(paste("theta_raw[", 1:10, "]"),
                  paste("theta[", 1:10, "]"),
                  "theta[100]", "lp__"), 
      probs=c(0.10, 0.5, 0.90));
