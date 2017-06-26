# fit WinBUGS model
library(R2OpenBUGS);
library(rstan)  
options(max.print = 5000);

source("carlin_data.R");

N = data$N;
inits <- function (N) {list(beta0 = rnorm(1, 0, 2),
                           beta1 = rnorm(1, 0, 2),
                           tau_theta = exp(runif(1, -2, 2)),
                           tau_phi = exp(runif(1, -2, 2)),
                           theta = rnorm(N, 0, 2),
                           phi = rep(0, N))}

inits=list(inits(N),inits(N),inits(N),inits(N));

params2save = c("beta0", "beta1", "tau_phi", "tau_theta", "sd_phi", "sd_theta", "psi", "eta[]", "phi[]", "theta[]");

pathOpenBugs = "/Users/mitzi/.wine/drive_c/Program\ Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe";

fit_bugs_200K <- bugs(data, inits, params2save,
     n.iter=200000,
     n.chains=4,
     n.burnin=190000,
     n.thin=2,
     model.file="carlin_bym.txt",
     digits=7,
     OpenBUGS.pgm=pathOpenBugs,
     useWINE=TRUE,
     WINE="/usr/local/bin/wine",
     debug = TRUE,
     bugs.seed=1);

capture.output(rstan::monitor(fit_bugs_200K$sims.array, digits=3,probs=c(0.025, 0.975), print=TRUE), file="scot_bugs_200K_iters.txt");

# fit_bugs_500K <- bugs(data, inits, params2save,
#      n.iter=500000,
#      n.chains=4,
#      n.burnin=490000,
#      n.thin=2,
#      model.file="carlin_bym.txt",
#      digits=7,
#      OpenBUGS.pgm=pathOpenBugs,
#      useWINE=TRUE,
#      WINE="/usr/local/bin/wine",
#      debug = TRUE,
#      bugs.seed=1);
# 
# capture.output(rstan::monitor(fit_bugs_500K$sims.array, digits=3,probs=c(0.025, 0.975), print=TRUE), file="scot_bugs_500K_iters.txt");
# 
# 
# fit_bugs_1M <- bugs(data, inits, params2save,
#      n.iter=1000000,
#      n.chains=4,
#      n.burnin=999000,
#      n.thin=2,
#      model.file="carlin_bym.txt",
#      digits=7,
#      OpenBUGS.pgm=pathOpenBugs,
#      useWINE=TRUE,
#      WINE="/usr/local/bin/wine",
#      debug = TRUE,
#      bugs.seed=1);
# 
# capture.output(rstan::monitor(fit_bugs_1M$sims.array, digits=3,probs=c(0.025, 0.975), print=TRUE), file="scot_bugs_1M_iters.txt");
# 
# 
# 
