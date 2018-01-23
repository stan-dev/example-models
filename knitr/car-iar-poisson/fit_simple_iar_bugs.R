library(R2OpenBUGS);
library(rstan)  

source("scotland_data.R");

iter = 100000;
burn =  90000;
mfile = "simple_iar.txt";
ofile = "simple_iar_bugs_100K_iters.txt";
thin = 2;
nchain = 2;

inits <- function () {list(phi = rep(0, 56))};
params2save = c("phi");

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
      debug=TRUE,
      WINE="/usr/local/bin/wine",
      DIC=FALSE,
      bugs.seed=1);

capture.output(rstan::monitor(fit_bugs$sims.array,
      digits=3,probs=c(0.025, 0.975), warmup=0, print=TRUE),
      file=ofile);

attach.bugs(fit_bugs);
