library(R2OpenBUGS);
library(rstan)  

source("scotland_data.R");


iter = 100000;
burn =  90000;
mfile = "bym_bugs.txt";
ofile = "bym_bugs_100K_iters.txt";
thin = 2;
nchain = 2;

data$x = 0.1 * data$x;
inits <- function () {
         list(beta0 = rnorm(1,0,1),
         beta1 = rnorm(1,0,1),
         tau_theta = exp(rnorm(1,0,1)),
         tau_phi = exp(rnorm(1,0,1)),
         theta = rnorm(56,0,1),
         phi = rep(0,56))};
 
params2save = c("beta0", "beta1", "sigma_phi", "sigma_theta",
  "tau_phi", "tau_theta", "mu", "phi", "theta");

pathOpenBugs = "/Users/mitzi/.wine/drive_c/Program\ Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe";

fit_bugs <- bugs(data, inits, params2save,
      n.iter=iter,
      n.chains=nchain,
      n.burnin=burn,
      n.thin=thin,
      model.file=mfile,
      digits=7,
      OpenBUGS.pgm=pathOpenBugs,
      useWINE=TRUE,
      DIC=FALSE,
      WINE="/usr/local/bin/wine",
      bugs.seed=1);

capture.output(rstan::monitor(fit_bugs$sims.array,
      digits=3,probs=c(0.025, 0.975), warmup=0, print=TRUE),
      file=ofile);

sims = rstan::monitor(fit_bugs$sims.array,
      probs=c(0.025, 0.975), warmup=0, print=TRUE);
