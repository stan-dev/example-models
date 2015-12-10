## 4. Introduction to random effects: Conventional Poisson GLMM for
## count data
## 4.2. Accounting for overdispersion by random effects-modeling in R
## and WinBUGS
## 4.2.1. Generation and analysis of simulated data

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Generate simualted data
## data.fn() is defined in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
set.seed(123)
data <- data.fn()
mny <- mean(data$year)
sdy <- sd(data$year)
cov1 <- (data$year - mny) / sdy

## Bundle data
stan_data <- list(C = data$C, n = length(data$C), year = cov1)

## Initial values
inits <- function() { list(alpha = runif(1, -2, 2),
                           beta1 = runif(1, -3, 3),
                           sd = runif(1, 0, 1))}

## Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda", "sigma",
            "eps")

# MCMC settings
ni <- 15000
nt <- 10
nb <- 5000
nc <- 4

## Call Stan from R
out <- stan("GLMM_Poisson.stan", data = stan_data,
            init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)

print(out)
