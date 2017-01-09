## 9. Multistate capture-recapture models
## 9.7. Real-data example: the showy lady's slipper

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data file "orchids.txt" is available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
rCH <- as.matrix(read.table("orchids.txt", sep = " ", header = FALSE))
rCH[rCH == 0] <- 3

stan_data <- list(y = rCH, n_occasions = dim(rCH)[2], nind = dim(rCH)[1])

## Parameters monitored
params <- c("s", "psiV", "psiF", "psiD")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(s = runif((dim(rCH)[2] - 1), 0, 1)))

## Call Stan from R
ls <- stan("ladyslipper.stan",
           data = stan_data, init = inits, pars = params,
           chains = nc, iter = ni, warmup = nb, thin = nt,
           seed = 1,
           open_progress = FALSE)
print(ls, digits = 3)
