## 12. Metapopulation modeling of abundance using
##     hierarchical Poisson regression: binomial mixture models
## 12.3. Analysis of real data: Open-population binomial-mixture model
## 12.3.3. Binomial-mixture model with overdispersion in both abundance
##         and detection


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(1234)

## Read data
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
params <- c("totalN", "alpha_lam", "beta", "sd_lam", "sd_p",
            "mean_abundance", "mean_N", "mean_detection",
            "fit", "fit_new")

## MCMC settings
ni <- 3500
nt <- 3
nb <- 500
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(alpha_lam = runif(7, -1, 1),
         sd_p = runif(1, 0, 3)))

## Call Stan from R
out2 <- stan("Nmix2.stan",
             data = list(y = y, R = R, T = T,
                         first = first, last = last, K = 70),
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 2,
             open_progress = FALSE)
print(out2, digits = 3)

## Note: This model will take a long time to complete. Some posteriors
## will differ from those in the book. This may be due to lower
## precision of them, or due to difference in prior of sd_lam to some
## extent. For simulated data, this model will produce similar
## results to BUGS.
