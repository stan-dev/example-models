library(rstan)
theme_set(theme_bw())

rstan_options(auto_write=TRUE)

set.seed(567583)

stan_pk_model <- stan_model("oral_1cmt_mm_run.stan")

source("oral_1cmt_mm_run.data.R", echo=TRUE)

ka <- log(2)/1
V <- 5
Km <- 3           ## conc at which 50% "saturation" occurs
k0 <- log(2)/8    ## elimination at very small conc
Vm <- V * Km * k0 ## zero order elimination when saturated
Am <- V * Km      ## amount at which 50% saturation occurs
theta <- log(c(ka, k0, Am, V))
theta_trans <- theta[1:3]
theta_trans[1:2] <- exp(theta_trans[1:2])

J <- max(id)

inits <- function() {
    list(theta_lelim=rnorm(2, theta[2:1], 0.5),
         theta_raw=rnorm(2, theta[3:4], 0.5),
         xi=matrix(theta[c(2,4)], 2, J),
         omega=rlnorm(2, log(0.1), 1),
         sigma_y=rlnorm(1, log(0.1), 1))
}


fit <- sampling(stan_pk_model, chains=1, warmup=200, iter=400, init=inits, refresh=1)

pars <- c("theta", "omega", "sigma_y")

print(fit, pars)
theta

traceplot(fit, pars)

pairs(fit, pars=pars)
