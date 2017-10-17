library(rstan)   
options(mc.cores = parallel::detectCores())  

source("mungeCARdata4stan.R")  
source("scotland_data.R")

iter = 10000;
mfile = "simple_iar.stan";
ofile = "simple_iar_stan_010K_iters.txt";

nbs = mungeCARdata4stan(data$adj, data$num);
N = data$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

fit_stan = stan(mfile, data=list(N,N_edges,node1,node2), iter=iter);

capture.output(print(fit_stan, digits=3, probs=c(0.025, 0.975)),
               file=ofile);

# estimate posterior co-variance for different parameters
#  neighbors:  phi[24]: 27,30,31,44,47,48,55,56,
#              phi[23]: 9,29,34,36,37,39,
#              phi[6]:  3, 8,
#              phi[8]:  6

stanfit = extract(fit_stan);
# cov neighbors
cov(stanfit$phi[,6],stanfit$phi[,8]);
cov(stanfit$phi[,6],stanfit$phi[,3]);
cov(stanfit$phi[,10],stanfit$phi[,22]);
# cov non-neighbors
cov(stanfit$phi[,6],stanfit$phi[,54]);
cov(stanfit$phi[,8],stanfit$phi[,54]);
cov(stanfit$phi[,2],stanfit$phi[,55]);
cov(stanfit$phi[,1],stanfit$phi[,55]);
cov(stanfit$phi[,2],stanfit$phi[,50]);

