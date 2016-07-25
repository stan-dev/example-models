library(rstan)
library(plyr)
theme_set(theme_bw())

set.seed(2345)

## 1cmt oral dosing pharmacokinetic model with multiple dosing (using
## Nonmem addl coding)

## simulation scenario: single loading dose (100% higher), then daily
## dosing for 6 days, 

## time unit is [h]

dose_lamt <- c(log(30), rep(log(15), 6))
dose_time <- 24 * seq(0, 6)
dose_cmt  <- rep(1, 7)
dose_tau  <- rep(0, 7)
dose_addl <- rep(0, 7)

init_lstate <- rep(-5, 2)
init_time <- 0

## log of ka, Vm, Km, V
ka <- log(2)/1
V <- 5
Km <- 3  ## conc at which 50% "saturation" occurs
k0 <- log(2)/8 ## elimination rate at very small conc
Vm <- V * Km * k0 ## maximal mass rate elimiationt at high conc
Am <- V * Km      ## mass at 50% of saturation
theta <- log(c(ka, k0, Am, V))
theta_trans <- theta[1:3]
theta_trans[1:2] <- exp(theta_trans[1:2])
lscale <- c(0, log(V))

## no of patients
J <- 20

## simulate per patient parameters of ke and V
Theta <- matrix(theta[1:3], J, 3, byrow=TRUE)
Lscale <- matrix(c(0, theta[4]), J, 2, byrow=TRUE)
Init_time <- rep(init_time, J)
Init_lstate <- matrix(-4, J, 2)

## simulate subject specific ke (at most 50% deviation)
omega_k0 <- log(1.5)/1.96

Theta[,2] <- rnorm(J, Theta[,2], omega_k0)

## ensure that we have no flip flop
all((Theta[,1] - Theta[,2]) > 0)
summary((Theta[,1] - Theta[,2]) )

## simulate subject specific V (at most 30% deviation)
omega_V <- log(1.3)/1.96

Lscale[,2] <- rnorm(J, Lscale[,2], omega_V)


Theta_trans <- Theta[,1:3]
Theta_trans[,1:2] <- exp(Theta_trans[,1:2])


## 10% relative residual error
sigma_y <- 0.1

## assemble Stan model
stan_pk_model <- stanc_builder("oral_1cmt_mm.stan")
## and export functions to R, used for simulation
expose_stan_functions(stan_pk_model)

## define design, i.e 6 observations on densely sampled days (first
## and last); in between we only measure the trough concentration
## (just before dosing)
Ndaily <- 6
obs_time_dense <- round(10^(seq(log10(0.1), log10(16), length=Ndaily)), 2)
obs_time_dense <- c(seq(0.5,7, length=Ndaily/2), seq(10, 15, length=Ndaily/2))
obs_time_dense

obs_time <- sort(c(seq(24, 24 * 6, by = 24),
                   obs_time_dense,
                   obs_time_dense + 24*6,
                   obs_time_dense + 24*8))

obs_time

x_r <- c(0)
x_i <- c(0)

## visualize population profile and design
curve(exp(pk_model(dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl, init_lstate, -1E-4, 24*x, theta_trans, lscale, x_r, x_i)[,2]), 0, 12, n=501, ylab="Conc", xlab="Time [d]")
points(obs_time/24, rep(1, length(obs_time)))

curve(pk_model(dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl, init_lstate, -1E-4, 24*x, theta_trans, lscale, x_r, x_i)[,2], 0.01, 12, n=501, ylab="Conc", xlab="Time [d]")
points(obs_time/24, rep(1, length(obs_time)))

## create Nonmem data set

## observation part for each patient
nm_obs <- data.frame(time=obs_time, cmt=2, evid=0, amt=0, tau=0, addl=0, mdv=0, dv=0)
## dosing part for each patient
nm_dose <- data.frame(time=dose_time+1E-5, cmt=1, evid=1, amt=exp(dose_lamt), tau=dose_tau, addl=dose_addl, mdv=1, dv=0)

nm <- arrange(rbind(nm_obs, nm_dose), time, evid, cmt)

## setup for all patients in the same way
copy_ind <- rep(1:nrow(nm), each=J)
nm <- nm[copy_ind,]
nm$id <- 1:J

## sort by id, time, evid (first observation, then dosing), finally cmt
nm <- arrange(nm, id, time, evid, cmt)

## simulate design
sim <- do.call(evaluate_model_nm, c(nm[,names(nm) != "dv"], list(init_lstate=Init_lstate, theta=Theta_trans, lscale=Lscale, init_time=Init_time-1E-3, x_r=x_r, x_i=x_i)))

## save simulated true mean conc in sdv
nm$sdv <- round(exp(sim[,2]), 3)
## apply noise
nm$dv <- round(rlnorm(nrow(nm), log(nm$sdv), sigma_y), 3)
nm$sdv[nm$cmt != 2] <- nm$dv[nm$cmt != 2] <- 0


qplot(time, dv, data=subset(nm, mdv==0), group=factor(id), geom="line", log="y")
qplot(time, dv, data=subset(nm, mdv==0), group=id, geom="line")

stan_data <- c(nm, list(N=nrow(nm),
                        prior_theta_mean=theta,
                        prior_theta_sd=log(c(10, 10, 10, 10))/1.96
                        )
               )

## save stan data
stan_rdump(names(stan_data), "oral_1cmt_mm_run.data.R", envir=list2env(stan_data))

## also serialize out Stan program with all includes applied
cat(stan_pk_model$model_code, file="oral_1cmt_mm_run.stan")

## finally save Nonmem data set
write.csv(nm, file="oral_1cmt_mm_run_nm.csv", quote=FALSE, row.names=FALSE)
