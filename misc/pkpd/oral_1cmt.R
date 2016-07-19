library(rstan)
theme_set(theme_bw())

rstan_options(auto_write=TRUE)

set.seed(567583)

stan_pk_model <- stan_model("oral_1cmt_run.stan")

source("oral_1cmt_run.data.R", echo=TRUE)

fit <- sampling(stan_pk_model, chains=2, warmup=500, iter=1000)

pars <- c("theta", "omega", "sigma_y")

print(fit, pars)

