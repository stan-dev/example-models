## 13. Metapopulation modeling of species distributions using
##     hierarchical logistic regression: Site-occupancy models
## 13.5. Dynamic (multi-season) site-occupancy models
## 13.5.2. Dynamic occupancy modeling in a real data set

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data file "burnet.txt" is available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
bdat <- read.table(file = "burnet.txt", header = T)
y <- array(NA, dim = c(95, 2, 7))   # 95 sites, 2 reps, 7 days
for (i in 1:7) {
  sel.rows <- bdat$day == i
  y[,,i] <- as.matrix(bdat)[sel.rows, 3:4]
}
y[y > 0] <- 1   # Convert counts to detection/nondetection data
# Aggregate detections over reps within a day and bundle data
yy <- apply(y, c(1, 3), sum, na.rm = TRUE)
stan_data <- list(y = yy, nsite = dim(yy)[1],
                  nyear = dim(yy)[2])

## Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n_occ", "growthr", "turnover")

## MCMC settings
ni <- 2000
nt <- 1
nb <- 1000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(psi1 = runif(1, 0, 1),
         p = runif(1, 0, 1)))

## Call Stan from R
out2 <- stan("Dynocc2.stan",
             data = stan_data,
             init = inits, pars = params,
             chains = nc, iter = ni, warmup = nb, thin = nt,
             seed = 1,
             control = list(adapt_delta = 0.8),
             open_progress = FALSE)
print(out2, digits = 3)

## Note: the original BUGS model in the BPA book is incorrect.
##       See errata in the website of the book.
