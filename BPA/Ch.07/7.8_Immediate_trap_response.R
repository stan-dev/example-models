## 7. Estimation of survival probabilities using capture-recapture data
## 7.8. Immediate trap response in recapture probability

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Import data
## The data file "trap.txt" is available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
CH <- as.matrix(read.table(file = "trap.txt", sep = " "))
# Create matrix m indicating when an individual was captured
m <- CH[, 1:(dim(CH)[2] - 1)]
u <- which(m == 0)
m[u] <- 2

## Bundle data
stan_data <- list(y = CH,
                  nind = dim(CH)[1],
                  n_occasions = dim(CH)[2],
                  m = m)

## Parameters monitored
params <- c("mean_phi", "beta")

## MCMC settings
ni <- 3000
nt <- 2
nb <- 1000
nc <- 4

## Initial values
inits <- function() list(mean_phi = runif(1, 0, 1),
                         beta = runif(2, 0, 1))

## Call Stan from R
cjs_trap  <- stan("cjs_trap.stan",
                  data = stan_data, init = inits, pars = params,
                  chains = nc, iter = ni, warmup = nb, thin = nt,
                  seed = 1,
                  open_progress = FALSE)

## Summarize posteriors
print(cjs_trap, digits = 3)
