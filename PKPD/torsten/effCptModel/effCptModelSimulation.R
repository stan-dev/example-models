## Template to simulate PKPD data
##
## Effect Compartment Model example from ACoP7 workshop

library(dplyr)
library(rstan)

modelName <- "effCptModel"

xdata1 <- read.csv("phase1effcpt.csv", as.is = TRUE)
xdata2 <- read.csv("phase2effcpt.csv", as.is = TRUE)

xdata1 <- xdata1 %>%
  select(id = subject, time, weight, dose, cobs, resp = fxa.inh.obs) %>%
  mutate(cobs = as.numeric(cobs),
         study = 1,
         evid = 0) %>%
  filter(dose > 0, time > 0)

dose1 <- xdata1 %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(time = 0,
         evid = 1,
         cobs = NA,
         resp = NA)

xdata1 <- xdata1 %>%
  bind_rows(dose1) %>%
  arrange(id, time, desc(evid))

xdata2 <- xdata2 %>%
  filter(drug == "ME-2") %>%
  select(id = patient, time, weight, dose, cobs, resp = fxa.inh.obs) %>%
  mutate(cobs = as.numeric(cobs),
         study = 2,
         evid = 0,
         id = id + max(xdata1$id)) %>%
  filter(time > 0)

dose2 <- xdata2 %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(evid = 1,
         cobs = NA,
         resp = NA) %>%
  select(-time) %>%
  merge(data.frame(time = seq(0, 13 * 12, by = 12)))

xdata2 <- xdata2 %>%
  bind_rows(dose2) %>%
  arrange(id, time, desc(evid))

xdata <- xdata1 %>%
  bind_rows(xdata2)

nt <- nrow(xdata)
start <- (1:nt)[!duplicated(xdata$id)]
end <- c(start[-1] - 1, nt)

## Indices of records containing observed concentrations
iObs <- with(xdata, (1:nrow(xdata))[!is.na(cobs) & evid == 0])
nObs <- length(iObs)

## Additional data required for LinCptModel
rate <- rep(0, nt)
ii <- rep(0, nt)
addl <- rep(0, nt)
ss <- rep(0, nt)

## create data set
data <- with(xdata,
             list(
               nSubjects = length(unique(id)),
               nt = nt,
               nObs = nObs,
               iObs = iObs,
               amt = dose,
               cmt = rep(1, nt),
               evid = evid,
               start = start,
               end = end,
               time = time,
               rate = rate,
               ii = ii,
               addl = addl,
               ss = ss,
               cObs = cobs[iObs],
               respObs = resp[iObs],
               weight = weight[!duplicated(id)]
             ))

## create initial estimates
init <- function(){
  list(CLHat = exp(rnorm(1, log(10), 0.2)),
       QHat = exp(rnorm(1, log(20), 0.2)),
       V1Hat = exp(rnorm(1, log(70), 0.2)),
       V2Hat = exp(rnorm(1, log(70), 0.2)),
       kaHat = exp(rnorm(1, log(2), 0.2)),
       ke0Hat = exp(rnorm(1,log(1),0.2)),
       EC50Hat = exp(rnorm(1,log(100),0.2)),
       omega = exp(rnorm(5, log(0.25), 0.5)),
       rho = diag(5),
       omegaKe0 = exp(rnorm(1, log(0.25), 0.5)),
       omegaEC50 = exp(rnorm(1, log(0.25), 0.5)),
       sigma = 0.5,
       sigmaResp = 20,
       logtheta = matrix(rep(log(c(10, 20, 70, 70, 2)), ea = 200), nrow = 200),
       logKe0 = rep(log(1), 200),
       logEC50 = rep(log(100), 200))
}

with(data, stan_rdump(ls(data), file = paste0(modelName, ".data.R")))
inits <- init()
with(inits, stan_rdump(ls(inits), file = paste0(modelName, ".init.R")))
