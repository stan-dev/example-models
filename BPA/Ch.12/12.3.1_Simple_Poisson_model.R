## 12. Metapopulation modeling of abundance using
##     hierarchical Poisson regression: binomial mixture models
## 12.3. Analysis of real data: Open-population binomial-mixture model
# 12.3.1. Simple Poisson model

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data file "fritilary.txt" is available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
bdat <- read.table("fritillary.txt", header = TRUE)
y <- array(NA, dim = c(95, 2, 7))	# 95 sites, 2 reps, 7 days

for(k in 1:7) {
  sel.rows <- bdat$day == k
  y[,,k] <- as.matrix(bdat)[sel.rows, 3:4]
}
R = nrow(y)
T = ncol(y)
first <- sapply(1:dim(y)[1], function(i)
    min(grep(FALSE, is.na(y[i, 1, ]))))
last <- sapply(1:dim(y)[1], function(i)
    max(grep(FALSE, is.na(y[i, 1, ]))))
y[is.na(y)] <- -1

## Parameters monitored
params <- c("totalN", "mean_abundance", "alpha_lam", "p", "fit", "fit_new")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(alpha_lam = runif(7, -1, 1)))

## Call Stan from R
out0 <- stan("Nmix0.stan",
             data = list(y = y, R = R, T = 2,
                         first = first, last = last, K = 100),
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             control = list(adapt_delta = 0.9),
             open_progress = FALSE)
print(out0, digits = 3)

## Note: Posteriors relevant to day 1 and 2 will somewhat differ from
## those in the book. This may be due to small number of observations
## in those days (total = 4 and 0, respectively), as described in the
## text. For simulated data, this model will produce similar results
## to BUGS.
