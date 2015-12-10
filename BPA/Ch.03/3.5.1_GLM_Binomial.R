## 3. Introduction to the generalized linear model (GLM): The simplest
## model for count data
## 3.5. Binomial GLM for modeling bounded counts or proportions
## 3.5.1. Generation and analysis of simulated data

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Generate simulated data
## data.fn() is defined in bpa-code.txt
set.seed(123)
data <- data.fn(nyears = 40, alpha = 1, beta1 = -0.03, beta2 = -0.9)

## Bundle data
stan_data <- list(C = data$C, N = data$N,
                  nyears = length(data$C), year = data$YR)

## Initial values
inits <- function() list(alpha = runif(1, -1, 1),
                         beta1 = runif(1, -1, 1),
                         beta2 = runif(1, -1, 1))

## Parameters monitored
params <- c("alpha", "beta1", "beta2", "p")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Call Stan from R
out <- stan("GLM_Binomial.stan", data = stan_data,
            init = inits, pars = params,
            chains = nc, thin = nt, iter = ni, warmup = nb,
            seed = 1,
            open_progress = FALSE)

print(out)
