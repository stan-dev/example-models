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


# estimate posterior co-variance for different parameters
#  neighbors:  phi[24]: 27,30,31,44,47,48,55,56,   
#              phi[23]: 9,29,34,36,37,39,
#              phi[6]:  3, 8,
#              phi[8]:  6

attach.bugs(fit_bugs);
# cov neighbors
cov(phi[,6],phi[,8]);
cov(phi[,6],phi[,3]);
cov(phi[,10],phi[,22]);
# cov non-neighbors
cov(phi[,6],phi[,54]);
cov(phi[,8],phi[,54]);
cov(phi[,2],phi[,55]);
cov(phi[,1],phi[,55]);
cov(phi[,2],phi[,50]);

