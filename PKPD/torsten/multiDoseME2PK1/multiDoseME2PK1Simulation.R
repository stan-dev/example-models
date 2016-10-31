## EXAMPLE FOR STAN TORSTEN
## See example 3 from MetrumBayesianStanPAGE2016

rm(list = ls())
gc()

modelName <- "multiDoseME2PK1"
library(rstan)
set.seed(11191951) ## not required but assures repeatable results

################################################################################################
## get data file
xdata <- read.csv("fxaNONMEMData.csv", as.is = TRUE)

xdata <- as.keyed(as.best(xdata), key = c("ID", "TIME", "CMT", "EVID"))

nt <- nrow(xdata)
start <- (1:nt)[!duplicated(xdata$ID)]
end <- c(start[-1] - 1, nt)

## Indices of records containing observed concentrations
iObs <- with(xdata, (1:nrow(xdata))[!is.na(DV) & EVID == 0])
nObs <- length(iObs)

## create data set
data <- with(xdata,
             list(
                 nSubjects = length(unique(ID)),
                 nt = nt,
                 nObs = nObs,
                 iObs = iObs,
                 amt = 1000 * AMT,
                 rate = RATE,
                 ii = II,
                 cmt = CMT,
                 evid = EVID,
                 addl = ADDL,
                 ss = SS,
                 start = start,
                 end = end,
                 time = TIME,
                 cObs = DV[iObs],
                 weight = WEIGHT[!duplicated(ID)]
))

## create initial estimates
init <- function(){
    list(CLHat = exp(rnorm(1, log(10), 0.2)),
         QHat = exp(rnorm(1, log(15), 0.2)),
         V1Hat = exp(rnorm(1, log(35), 0.2)),
         V2Hat = exp(rnorm(1, log(105), 0.2)),
         kaHat = exp(rnorm(1, log(2), 0.2)),
         omega = exp(rnorm(5, log(0.25), 0.5)),
         rho = diag(5),
         sigma = runif(1, 0.5, 2),
         logtheta = matrix(rep(log(c(10, 15, 35, 105, 2)), ea = 212), nrow = 212))
}

with(data, stan_rdump(ls(data), file = paste0(modelName, ".data.R")))
inits <- init()
with(inits, stan_rdump(ls(inits), file = paste0(modelName, ".init.R")))
