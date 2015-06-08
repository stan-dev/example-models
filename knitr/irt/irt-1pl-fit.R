# ---- fit-1pl-stan-compile ----
library(rstan);
model <- stan_model("irt_1pl.stan");

## ---- fit-1pl-stan-init ----
init_fun <- function(chain_id) {
  return(list(theta=runif(J, -2, 2), b=runif(I, -2, 2)));
}

## ---- fit-1pl-stan-sampling ----
fit <-sampling(model, data = c("I", "J", "y"), init=init_fun, refresh=2000, seed=1234)

## ---- fit-1pl-stan-print ----
options("width"=100);
print(fit, c(paste("b[", 1:10, "]"), paste("theta[", 1:10, "]"), "lp__"), 
      probs=c(0.1, 0.5, 0.9));

## ---- fit-1pl-stan-compare ----
print(b[1:10], digits=2);
print(theta[1:10], digits=2);

