## ---- fit-1pl-stan-pin ----
model_pin <- stan_model("irt_1pl_pin.stan");
fit_pin <- sampling(model_pin, data=c("I", "J", "y"), refresh=2000, seed=1234);

## ---- fit-print-1pl-stan-pin ----
print(fit_pin, c(paste("b[", 1:10, "]"), paste("theta[", 1:10, "]"), "lp__"), 
      probs=c(0.10, 0.5, 0.90));
print("theta[1:10] + theta[100] =", quote=FALSE)
print(theta[1:10] + theta[100], quote=FALSE, digits=1);
