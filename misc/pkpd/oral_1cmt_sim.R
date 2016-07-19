library(rstan)
library(plyr)
theme_set(theme_bw())

set.seed(2345)

## 1cmt oral dosing pharmacokinetic model with multiple dosing (using
## Nonmem addl coding)

## simulation scenario: single loading dose (50% higher), then daily
## dosing for 6 days, 

## time unit is [h]

dose_lamt <- c(log(15), log(10))
dose_time <- c(0, 24)
dose_cmt  <- c(1, 1)
dose_tau  <- c(0, 24)
dose_addl <- c(0, 5) ## these are 5 additional doses (in total 6)

init_lstate <- rep(-25, 2)
init_time <- 0

## log of ka, ke, V
theta <- log(c(log(2)/2, log(2)/12, 10))
lscale <- rep(0, 2)

## no of patients
J <- 30

## simulate per patient parameters of ke and V
Theta <- matrix(theta[1:2], J, 2, byrow=TRUE)
Lscale <- matrix(c(0, theta[3]), J, 2, byrow=TRUE)
Init_time <- rep(0, J)
Init_lstate <- matrix(-25, J, 2)

## simulate subject specific ke (at most 50% deviation)
omega_ke <- log(1.5)/1.96

Theta[,2] <- rnorm(J, Theta[,2], omega_ke)

## simulate subject specific V (at most 30% deviation)
omega_V <- log(1.3)/1.96

Lscale[,2] <- rnorm(J, Lscale[,2], omega_V)

## 10% relative residual error
sigma_y <- 0.1

## assemble Stan model
stan_pk_model <- stanc_builder("oral_1cmt.stan")
## and export functions to R, used for simulation
expose_stan_functions(stan_pk_model)

## define design, i.e 6 observations on densely sampled days (first
## and last); in between we only measure the trough concentration
## (just before dosing)
Ndaily <- 6
obs_time_dense <- round(10^(seq(log10(0.1), log10(16), length=Ndaily)), 2)
obs_time_dense <- c(seq(0.5,7, length=Ndaily/2), seq(10, 15, length=Ndaily/2))
obs_time_dense

obs_time <- sort(c(seq(24, 24 * 6, by = 24), obs_time_dense, obs_time_dense + 24*6))

obs_time

## visualize population profile and design
curve(exp(pk_model(dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl, init_lstate, 0, 24*x, theta, lscale)[,2]), 0, 8, n=501, ylab="Conc", xlab="Time [d]")
points(obs_time/24, rep(1, length(obs_time)))

## create Nonmem data set

## observation part for each patient
nm_obs <- data.frame(time=obs_time, cmt=2, evid=0, amt=0, tau=0, addl=0, mdv=0, dv=0)
## dosing part for each patient
nm_dose <- data.frame(time=dose_time, cmt=1, evid=1, amt=exp(dose_lamt), tau=dose_tau, addl=dose_addl, mdv=1, dv=0)

nm <- arrange(rbind(nm_obs, nm_dose), time, evid, cmt)

## setup for all patients in the same way
copy_ind <- rep(1:nrow(nm), each=J)
nm <- nm[copy_ind,]
nm$id <- 1:J

## sort by id, time, evid (first observation, then dosing), finally cmt
nm <- arrange(nm, id, time, evid, cmt)

## simulate design
sim <- do.call(evaluate_model_nm, c(nm[,names(nm) != "dv"], list(init_lstate=Init_lstate, theta=Theta, lscale=Lscale, init_time=Init_time)))

## save simulated true mean conc in sdv
nm$sdv <- round(exp(sim[,2]), 3)
## apply noise
nm$dv <- round(rlnorm(nrow(nm), log(nm$sdv), sigma_y), 3)
nm$sdv[nm$cmt != 2] <- nm$dv[nm$cmt != 2] <- 0


qplot(time, dv, data=subset(nm, mdv==0), group=factor(id), geom="line", log="y")
qplot(time, dv, data=subset(nm, mdv==0), group=id, geom="line")

stan_data <- c(nm, list(N=nrow(nm),
                        prior_theta_mean=log(c(log(2)/1, log(2)/10, 10)),
                        prior_theta_sd=log(c(10, 10, 10))/1.96
                        )
               )

## save stan data
stan_rdump(names(stan_data), "oral_1cmt_run.data.R", envir=list2env(stan_data))

## also serialize out Stan program with all includes applied
cat(stan_pk_model$model_code, file="oral_1cmt_run.stan")

## finally save Nonmem data set
write.csv(nm, file="oral_1cmt_run_nm.csv", quote=FALSE, row.names=FALSE)
