# fit WinBUGS model
library(R2OpenBUGS);

source("carlin_data.R");
source("carlin_inits.R");
source("carlin_inits_1.R");
source("carlin_inits_2.R");

inits=list(init_params,init_params_1,init_params_2);
params2save = c("beta0", "beta1", "tau.h", "tau.c", "sd.h", "sd.c", "alpha",
"theta[]","phi[]","xi[]");
pathOpenBugs = "/Users/mitzi/.wine/drive_c/Program\ Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe";

fit_bugs <- bugs(data, inits, params2save,
     n.iter=200000,
     n.chains=3,
     n.burnin=190000,
     model.file="carlin_bym.txt",
     digits=7,
     OpenBUGS.pgm=pathOpenBugs,
     useWINE=TRUE,
     WINE="/usr/local/bin/wine",
     bugs.seed=1);

library(rstan)  
print(rstan::monitor(fit_bugs$sims.array), digits=4) 
