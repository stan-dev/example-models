## 7. Estimation of survival probabilities using capture-recapture data
## 7.9. Parameter identifiability

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("cjs_t_t.data.R")

## Parameters monitored
params <- c("phi_t", "p_t")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i) {
    list(phi_t = runif((dim(stan_data$y)[2] - 1), 0, 1),
         p_t = runif((dim(stan_data$y)[2] - 1), 0, 1))})

## Call Stan from R
cjs_t_t  <- stan("cjs_t_t.stan",
                 data = stan_data, init = inits, pars = params,
                 chains = nc, iter = ni, warmup = nb, thin = nt,
                 seed = 1,
                 open_progress = FALSE)

## Summarize posteriors
print(cjs_t_t, digits = 3)
