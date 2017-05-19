## Template to simulate PKPD data
## Fribgerg-Karlsson population model

modelName <- "fribergKarlsson"

library(mrgsolve)
library(rstan)

## Simulate ME-2 plasma concentrations and ANC values
## using mrgsolve.

nSub <- 15; # number of subjects
nIIV <- 7; # number of parameters with inter-individual variations

code <- '
$PARAM CL = 10, Q = 15, VC = 35, VP = 105, KA = 2.0, MTT = 125, 
Circ0 = 5, alpha = 3E-4, gamma = 0.17, WT = 70

$SET delta=0.1 // simulation grid

$CMT GUT CENT PERI PROL TRANSIT1 TRANSIT2 TRANSIT3 CIRC

$MAIN
// Individual PK parameters
double CLi = exp(log(CL) + 0.75*log(WT/70) + ETA(1));
double Qi = exp(log(Q) + 0.75*log(WT/70) + ETA(2));
double VCi = exp(log(VC) + log(WT/70) + ETA(3));
double VPi = exp(log(VP) + log(WT/70) + ETA(4));

// Indivudal PD parameters
double MTTi = exp(log(MTT) + ETA(5));
double Circ0i = exp(log(Circ0) + ETA(6));
double alphai = exp(log(alpha) + ETA(7));

// Reparametrization
double k10 = CLi / VCi;
double k12 = Qi / VCi;
double k21 = Qi / VPi;
double ktr = 4/MTTi;

$ODE 
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - (k10 + k12) * CENT + k21 * PERI;
dxdt_PERI = k12 * CENT - k21 * PERI;
dxdt_PROL = ktr * (PROL + Circ0i) * 
((1 - alphai * CENT/VCi) * pow(Circ0i/(CIRC + Circ0i),gamma) - 1);
dxdt_TRANSIT1 = ktr * (PROL - TRANSIT1);
dxdt_TRANSIT2 = ktr * (TRANSIT1 - TRANSIT2);
dxdt_TRANSIT3 = ktr * (TRANSIT2 - TRANSIT3);
dxdt_CIRC = ktr * (TRANSIT3 - CIRC);

$OMEGA name="IIV"
0.0625 0.0625 0.0625 0.0625 0.04 0.0256 0.0256

$SIGMA 0.01 0.01 

$TABLE
double CP = CENT/VCi;
double DV1 = CENT/VCi * exp(EPS(1));
double DV2 = (CIRC + Circ0i) * exp(EPS(2));
double WEIGHT = WT;

$CAPTURE CP DV1 DV2 WEIGHT
'

mod <- mread("acum", tempdir(), code)
e1 <- expand.ev(amt = 80 * 1000, addl = 14, ii = 12, WT = rnorm(nSub, 70, 15))
out <- mod %>% data_set(e1) %>% carry.out(dose) %>% Req(CP,DV1,DV2) %>% mrgsim(end=500)

## Observation and dosing times
# doseTimes <- seq(0, 168, by = 12)
xpk <- c(0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8)
xpk <- c(xpk, xpk + 12, seq(24, 156, by = 12), c(xpk, 12, 18, 24) + 168)
xneut <- seq(0, 672, by = 48)
time <- sort(unique(c(xpk, xneut)))

## Assemble data set for Stan
xdata <- mod %>% data_set(e1) %>%
  carry.out(cmt, ii, addl, rate, amt, evid, ss) %>%
  mrgsim(Req = "DV1, DV2, WEIGHT", end = -1, add = time, rescort = 3) %>%
  as.data.frame

xdata <- xdata %>%
  mutate(DV1 = ifelse(time %in% xpk & time != 0 & evid == 0, DV1, NA),
         DV2 = ifelse(time %in% xneut & evid == 0, DV2, NA))

xdata$cmt[xdata$cmt == 0] <- 2 # switch from mrgsolve to NONMEM convention

nt <- nrow(xdata)

## Subject specific data
start <- (1:nt)[!duplicated(xdata$ID)]
end <- c(start[-1] -1, nt)
xsub <- subset(xdata, !duplicated(ID))
weight <- xsub$WEIGHT

## Look at simulated data using plots
p1 <- ggplot(xdata %>% filter(!is.na(DV1)), aes(x = time, y = DV1))
p1 + geom_point() + geom_line() +
  labs(x = "time (h)", y = "ME-2 plasma concentration (ng/mL)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) +
  facet_wrap(~ID)

p1 <- ggplot(xdata %>% filter(!is.na(DV2)), aes(x = time, y = DV2))
p1 + geom_point() + geom_line() +
  labs(x = "time (h)",
       y = "ANC") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) +
  facet_wrap(~ID)

## Indices of records containing observed concentrations
iObsPK <- with(xdata, (1:nrow(xdata))[!is.na(DV1) & evid == 0])
nObsPK <- length(iObsPK)
## Indices of records containing observed neutrophil counts
iObsPD <- with(xdata, (1:nrow(xdata))[!is.na(DV2) & evid == 0])
nObsPD <- length(iObsPD)

## Parameters for priors (only used for informed initial estimates)
CLHatPrior = 10
QHatPrior = 15
V1HatPrior = 35
V2HatPrior = 105
kaHatPrior = 2
CLHatPriorCV = 0.10
QHatPriorCV = 0.18
V1HatPriorCV = 0.14
V2HatPriorCV = 0.17
kaHatPriorCV = 0.16

circ0HatPrior <- 5
circ0HatPriorCV <- 0.20
mttHatPrior <- 125
mttHatPriorCV <- 0.2
gammaPrior <- 0.17
gammaPriorCV <- 0.2
alphaHatPrior <- 2.0E-4
alphaHatPriorCV <- 0.2

## create data set
data <- with(xdata,
             list(
               nt = nt,
               nObsPK = nObsPK,
               iObsPK = iObsPK,
               nObsPD = nObsPD,
               iObsPD = iObsPD,
               amt = amt,
               cmt = cmt,
               evid = evid,
               time = time,
               ii = ii,
               addl = addl,
               ss = ss,
               rate = rate,
               cObs = DV1[iObsPK],
               neutObs = DV2[iObsPD],
               
               nSubjects = nSub,
               nIIV = nIIV,
               start = start,
               end = end,
               weight = weight,
               
               # Priors for PD parameters
               circ0HatPrior = circ0HatPrior,
               circ0HatPriorCV = circ0HatPriorCV,
               mttHatPrior = mttHatPrior,
               mttHatPriorCV = mttHatPriorCV,
               gammaPrior = gammaPrior,
               gammaPriorCV = gammaPriorCV,
               alphaHatPrior = alphaHatPrior,
               alphaHatPriorCV = alphaHatPriorCV
             ))

## create initial estimates
init <- function(){
  list(CLHat = abs(rnorm(1, 0, 20)),
       QHat = abs(rnorm(1, 0, 20)),
       V1Hat = abs(rnorm(1, 0, 100)),
       V2Hat = abs(rnorm(1, 0, 1000)),
       kaHat = abs(rnorm(1, 0, 5)),
       sigma = 0.2,
       alphaHat = exp(rnorm(1, log(alphaHatPrior), alphaHatPriorCV)),
       mttHat = exp(rnorm(1, log(mttHatPrior), mttHatPriorCV)),
       circ0Hat = exp(rnorm(1, log(circ0HatPrior), circ0HatPriorCV)),
       gamma = exp(rnorm(1, log(gammaPrior), gammaPriorCV)),
       sigmaNeut = 0.2,
       L = diag(nIIV),
       omega = exp(rnorm(nIIV, log(0.05), 0.5)),
       etaStd = matrix(rep(0, nIIV * nSub), nrow = nIIV))
}

with(data, stan_rdump(ls(data), file = paste0(modelName,".data.R")))
inits <- init()
with(inits, stan_rdump(ls(inits), file = paste0(modelName,".init.R")))
