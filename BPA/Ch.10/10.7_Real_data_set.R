## 10. Estimation of survival, recruitment and population size using
##     the Jolly-Seber model
## 10.7. Analysis of a real data set: survival, recruitment and
##       population size of Leisler’s bats
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data file "leisleri.txt" is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
leis <- as.matrix(read.table("leisleri.txt", sep = " ", header = FALSE))
nz <- 300
CH.aug <- rbind(leis, matrix(0, ncol = dim(leis)[2], nrow = nz))
stan_data <- list(y = CH.aug, n_occasions = dim(CH.aug)[2],
                  M = dim(CH.aug)[1])

## Parameters monitored
params <- c("psi", "mean_p", "sigma2", "mean_phi", "N", "Nsuper", "b", "B")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(mean_phi = runif(1, 0, 1),
         mean_p = runif(1, 0, 1),
         sigma = runif(1, 0, 1)))

## Call Stan from R
nl <- stan("js_tempran.stan",
           data = stan_data, init = inits, pars = params,
           chains = nc, iter = ni, warmup = nb, thin = nt,
           seed = 1,
           open_progress = FALSE)
print(nl, digits = 3)
