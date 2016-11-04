## Template to simulate PKPD data

rm(list = ls())
gc()

modelName <- "LinTwoCptModel"

library(rstan)
library(mrgsolve)

set.seed(11091989) ## not required but assures repeatable results

###############################################################################################
## Simulate data using mrgsolve

code <- '
$PARAM CL=5, Q=8, V2=20, V3=70, KA=1.2

$CMT GUT CENT PERI

$GLOBAL
#define CP (CENT/V2)

$PKMODEL ncmt = 2, depot = TRUE

$SIGMA 0.025

$MAIN
pred_CL = CL;
pred_Q = Q;
pred_VC = V2;
pred_VP = V3;
pred_KA = KA;

$CAPTURE CP

$TABLE table(DV) = CP*exp(EPS(1));
'

mod <- mread("accum", tempdir(),code) %>% Req(DV) %>% update(end=480,delta=0.1)

e1 <- ev(amt=5000, rate = 150) # Create an initial dosing event
mod %>% ev(e1) %>% mrgsim(end=20) %>% plot # plot data

# create time at which data will be observed 
# NOTE: end time at t=20 -- if we go further, we get to the limit where CP -> 0. 
t1 <- seq(2,20,2)
t2 <- seq(0.25,2,0.25)
tall <- sort(c(t1,t2))

# save data in data frame 
SimData <- 
  mod %>% 
  ev(e1) %>% 
  carry.out(cmt,ii,addl,rate,amt,evid,ss) %>%
  mrgsim(Req="CP", end=-1, add=tall,recsort=3) %>%
  as.data.frame

SimData$cmt[SimData$cmt == 0] <- 2 ## adjust cmt (adopt NONMEM convention)
SimData <- SimData[!((SimData$evid == 0)&(SimData$DV == 0)),] ## remove observation with 0 drug concentration

################################################################################################
# Format Data for Stan using RStan

nt <- nrow(SimData)
iObs <- with(SimData, (1:nrow(SimData))[evid == 0])
nObs <- length(iObs)

## create Stan data set
data <- with(SimData,
             list(nt = nt,
                  nObs = nObs,
                  iObs = iObs,
                  time = time,
                  cObs = DV[iObs],
                  amt =  amt,
                  rate = rate,
                  cmt = cmt,
                  evid = evid,
                  ii = ii,
                  addl = addl,
                  ss = ss))

## create initial estimates
init <- function() 
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(20), 0.2)),
       V1 = exp(rnorm(1, log(70), 0.2)),
       V2 = exp(rnorm(1, log(70), 0.2)),
       ka = exp(rnorm(1, log(1), 0.2)),
       sigma = runif(1, 0.5, 2))

with(data, stan_rdump(ls(data), file = paste0(modelName, ".data.R")))
inits <- init()
with(inits, stan_rdump(ls(inits), file = paste0(modelName, ".init.R")))
