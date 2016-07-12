## EXAMPLE FOR STAN TORSTEN 
##
## Model Parameters: CL=5, Q=8, V1=20, V2=70, KA=1.2, SIGMA = 0.0025
## Single Dose 
##
## NOTE: can pick between TwoCptModel and GenTwoCptModel

rm(list = ls())
gc()

## Pick a model
modelName <- "TwoCptModelExample" # root names of model file
modelName <- "GenTwoCptModelExample"

diagnose <- FALSE # Set this to TRUE to run a diagnose instead of fitting a model
                  # A diagnose computes the adjoint of each parameter with respect to
                  # the posterior distribution using both automatic and finite differentiation. 

## Relative paths assuming the working directory is the script directory
## containing this script
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, "data", "derived")
modelDir <- file.path(projectDir, "model")
toolsDir <- file.path(scriptDir, "tools")
stanDir <- file.path(scriptDir, "cmd_lib", "cmdstan")
tempDir <- file.path(scriptDir, "temp")

dir.create(tempDir)

source(file.path(scriptDir, "pkgSetup.R"))
# source("installMrgsolve.R")

library(rstan)
library(parallel)
library(plyr)
library(dplyr)
library(metrumrg)
library(ggplot2)
library(mrgsolve)

source(file.path(toolsDir, "stanTools.R"))
source(file.path(toolsDir, "cmdStanTools.R"))

set.seed(11091989) ## not required but assures repeatable results

###############################################################################################
## Simulate data using mrgsolve

code <- '
$PARAM CL=5, Q=8, VC=20, VP=70, KA=1.2

$CMT GUT CENT PERI

$GLOBAL
#define CP (CENT/VC)

$ADVAN4 // Two compartment model

$SIGMA 0.0025

$MAIN
pred_CL = CL;
pred_Q = Q;
pred_VC = VC;
pred_VP = VP;
pred_KA = KA;

$CAPTURE CP

$TABLE table(DV) = CP*exp(EPS(1));
'

mod <- mread("accum", tempdir(),code) %>% Req(GUT, CENT, CP, DV) %>% update(end=480,delta=0.1)

e1 <- ev(amt=5000) # Create an initial dosing event
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
SimData <- SimData[!((SimData$evid == 0)&(SimData$CP == 0)),] ## remove observation with 0 drug concentration

################################################################################################

xdata <- SimData 

nt <- nrow(xdata)

iObs <- with(xdata, (1:nrow(xdata))[evid == 0])
nObs <- length(iObs)

## create Stan data set
data <- with(xdata,
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


## Specify the variables for which you want history and density plots
parametersToPlot <- c("CL", "Q", "V1", "V2", "ka", "sigma")

## Additional variables to monitor
otherRVs <- c("cObsPred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

################################################################################################
# run Stan
n.chains <- 4 # 4
nPost <- 1000 # 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1000 # 1000 ## Number of burn-in samples per chain after thinning
n.thin <- 1

n.iter <- nPost * n.thin
n.burnin <- nBurn * n.thin

with(data, stan_rdump(ls(data), file = file.path(tempDir, "data.R")))
inits <- init()
with(inits, stan_rdump(ls(inits), file = file.path(tempDir, "init.R")))

compileModel(model = file.path(modelDir, modelName), stanDir = stanDir)

RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()

chains <- 1:n.chains

if(diagnose == TRUE){
  ## Run diagnose: compare auto diff. to finite diff.
  mclapply(chains,
           function(chain, model, data, init)
             runDiagnose(model = model, data = data,
                         init = init, seed = sample(1:99999, 1),
                         chain = chain,
                         refresh = 100),
           model = file.path(modelDir, modelName),
           data = file.path(tempDir, "data.R"),
           init = file.path(tempDir, "init.R"),
           mc.cores = min(n.chains, detectCores()))
}

if(diagnose == FALSE){
  mclapply(chains,
           function(chain, model, data, iter, warmup, thin, init)
             runModel(model = model, data = data,
                      iter = iter, warmup = warmup, thin = thin,
                      init = init, seed = sample(1:999999, 1),
                      chain = chain,
                      ##                      max_depth = 10,
                      ##                      adapt_delta = 0.8,
                      ##                      stepsize = 1,
                      refresh = 100),
           model = file.path(modelDir, modelName),
           data = file.path(tempDir, "data.R"),
           init = file.path(tempDir, "init.R"),
           iter = n.iter, warmup = n.burnin, thin = n.thin,
           mc.cores = min(n.chains, detectCores()))
  
  
  posterior <- ldply(chains,
                     function(chain, modelName, modelDir){
                       res <- read.csv(paste(file.path(modelDir, modelName, modelName), chain, ".csv", sep = ""),
                                       as.is = TRUE, comment.char = "#")
                       cbind(chain = rep(chain, nrow(res)), iteration = 1:nrow(res), res)
                     }, modelName = modelName, modelDir = modelDir)
  
  ################################################################################################
  ## posterior distributions of parameters
  
  dir.create(figDir)
  dir.create(tabDir)
  
  graphics.off()
  ## open graphics device
  pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
      width = 6, height = 6, onefile = F)
  
  parametersToPlot <- names(posterior)[unlist(sapply(c(paste("^",parametersToPlot,"$",sep=""),
                                                       paste("^",parametersToPlot,".",sep="")),
                                                     grep, x = names(posterior)))]
  
  mcmcHistory(subset(posterior, select = c("chain", "iteration", parametersToPlot))) 
  pairs(subset(posterior, select = c("chain", "iteration", parametersToPlot)))
  
  mcmcDensity(subset(posterior, select = c("chain", "iteration", parametersToPlot)), byChain = TRUE)
  mcmcDensity(subset(posterior, select = c("chain", "iteration", parametersToPlot)))
  
  ptable <- parameterTable(subset(posterior, select = c("chain", "iteration", parametersToPlot)))
  write.csv(ptable, file.path(tabDir, paste(modelName, ".summary.csv", sep = "")))
  
  
  ################################################################################################
  ## posterior predictive distributions of plasma concentrations
  
  ## prediction of future observations in the same subjects, i.e., posterior predictions
  ## conditioned on observed data from the same subject
  
  pred <- as.matrix(posterior[, grep("cObsPred", dimnames(posterior)[[2]])])
  
  predSummary <- t(apply(pred,2,quantile,probs=c(0.05,0.5,0.95), na.rm = TRUE))
  colnames(predSummary) <- c("y5", "y50", "y95")
  predSummary <- xdata %>% filter(evid == 0) %>% bind_cols(as.data.frame(predSummary))
  
  p1 <- ggplot(predSummary, aes(x = time, y = DV))
  p1 <- p1 + geom_point() +
    labs(title = "individual predictions", x = "time (h)",
         y = "concentration") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none")
  print(p1 + geom_line(aes(x = time, y = y50)) +
          geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.25))
  
  dev.off()
}

