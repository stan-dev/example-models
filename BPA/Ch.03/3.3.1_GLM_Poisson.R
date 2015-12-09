## 3. Introduction to the generalized linear model (GLM): The simplest
## model for count data
## 3.3. Poisson GLM in R and WinBUGS for modeling times series of
## counts
## 3.3.1. Generation and analysis of simulated data

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Generate simulated data
## data.fn() is defined in bpa-code.txt
set.seed(123)
data <- data.fn()

## Bundle data
mean.year <- mean(data$year)             # Mean of year covariate
sd.year <- sd(data$year)                 # SD of year covariate

stan_data <- list(C = data$C, n = length(data$C),
                  year = (data$year - mean.year) / sd.year)

## Initial values
inits <- function() list(alpha = runif(1, -2, 2),
                         beta1 = runif(1, -3, 3))

## Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Stan
out <- stan("GLM_Poisson.stan", data = stan_data,
            init = inits, pars = params,
            chains = nc, thin = nt, iter = ni, warmup = nb,
            seed = 1,
            open_progress = FALSE)

print(out)
