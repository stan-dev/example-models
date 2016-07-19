library(rstan)
theme_set(theme_bw())

rstan_options(auto_write=TRUE)

set.seed(567583)

stan_pk_model <- stan_model("oral_1cmt_run.stan")

source("oral_1cmt_run.data.R", echo=TRUE)

theta <- log(c(log(2)/2, log(2)/12, 10))
J <- max(id)

inits <- function() {
    list(theta=rnorm(3, theta, 0.5),
         xi=matrix(theta[2:3], 2, J),
         omega=rlnorm(2, log(0.1), 1),
         sigma_y=rlnorm(1, log(0.1), 1))
}


fit <- sampling(stan_pk_model, chains=4, warmup=500, iter=1000, init=inits, refresh=10)

pars <- c("theta", "omega", "sigma_y")

print(fit, pars)

