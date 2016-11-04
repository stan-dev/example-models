## Template to simulate PKPD data for population model
## Inter-individual variability on all PK parammeters

rm(list = ls())
gc()

modelName <- "TwoCptModelPopulation"

library(rstan)
library(mrgsolve)

###############################################################################################
## mrgsolve

nSub <- 10; # number of subjects
nIIV <- 5; # number of parameters with Inter-Individual variation

code <- '
$PARAM CL=5, Q=8, V2=20, V3=70, KA=1.2

$SET delta=0.1 // simulation grid

$CMT GUT CENT PERI

$PKMODEL ncmt = 2, depot = TRUE

$MAIN
double CLi = exp(log(CL) + ETA(1));
double Qi = exp(log(Q) + ETA(2));
double V2i = exp(log(V2) + ETA(3));
double V3i = exp(log(V3) + ETA(4));
double KAi = exp(log(KA) + ETA(5));

pred_CL = CLi;
pred_Q = Qi;
pred_V2 = V2i;
pred_V3 = V3i;
pred_KA = KAi;

$OMEGA name="IIV"
0.0025 0.0025 0.0025 0.0025 0.0025

$SIGMA 0.025

$TABLE
table(DV) = CENT/V2*exp(EPS(1));
'

mod <- mread("accum", tempdir(),code)

e1 <- expand.ev(amt=rep(1000, nSub)) # Create an initial dosing event
out <- mod %>% data_set(e1) %>% carry.out(dose) %>% Req(DV) %>% mrgsim(end=50)
plot(out, DV~time|factor(ID),scales="same")

# create time at which data will be observed 
t1 <- seq(2,20,2)
t2 <- seq(0.25,2,0.25)
tall <- sort(c(t1,t2))

# save data in data frame 
SimData <- 
  mod %>%
  data_set(e1) %>%
  carry.out(cmt,ii,addl,rate,amt,evid,ss) %>%
  mrgsim(Req="DV", end=-1, add=tall, recsort=3) %>%
  as.data.frame

SimData$cmt[SimData$cmt == 0] <- 2 ## adjust cmt (adopt NONMEM convention)
SimData <- SimData[!((SimData$evid == 0)&(SimData$DV == 0)),] ## remove observation with 0 drug concentration

################################################################################################
## Format data for Stan 

nt <- nrow(SimData)

iObs <- with(SimData, (1:nrow(SimData))[evid == 0])
nObs <- length(iObs)

## Subject specific data
xsub <- subset(SimData, !duplicated(ID))
nSubjects <- length(xsub$ID)

## Row indices for start and end of each individual's data
start <- (1:nt)[!duplicated(SimData$ID)]
end <- c(start[-1] - 1, nt)

## create Stan data set
data <- with(SimData,
             list(nt = nt,
                  nSubjects = nSubjects,
                  start = start,
                  end = end,
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
                  ss = ss,
                  nIIV = nIIV)) ## number of parameters with IIV

## create initial estimates
init <- function()
  list(CLHat = exp(rnorm(1, log(10), 0.2)),
       QHat = exp(rnorm(1, log(20), 0.2)),
       V1Hat = exp(rnorm(1, log(70), 0.2)),
       V2Hat = exp(rnorm(1, log(70), 0.2)),
       kaHat = exp(rnorm(1, log(1), 0.2)),
       sigma = runif(1, 0.5, 2),
       L = diag(nIIV),
       etaStd = matrix(rep(0, nIIV * data$nSubjects), nrow = nIIV),
       omega = runif(nIIV, 0.5, 2),
       logtheta = matrix(rep(log(c(exp(rnorm(1, log(10), 0.2)),
                                   exp(rnorm(1, log(20), 0.2)),
                                   exp(rnorm(1, log(70), 0.2)),
                                   exp(rnorm(1, log(70), 0.2)),
                                   exp(rnorm(1, log(1), 0.2)))),
                             ea = data$nSubjects),
                         nrow = data$nSubjects))

with(data, stan_rdump(ls(data), file = paste0(modelName,".data.R")))
inits <- init()
with(inits, stan_rdump(ls(inits), file = paste0(modelName,".init.R")))
