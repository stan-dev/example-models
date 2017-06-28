library(rstan)   
options(mc.cores = parallel::detectCores())  

source("mungeCARdata4stan.R")  
source("carlin_data.R")

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

flist = extract(car_stanfit);
# cov neighbors
cov(flist$phi[,6],flist$phi[,8]);
cov(flist$phi[,6],flist$phi[,3]);
cov(flist$phi[,10],flist$phi[,22]);
# cov non-neighbors
cov(flist$phi[,6],flist$phi[,54]);
cov(flist$phi[,8],flist$phi[,54]);
cov(flist$phi[,2],flist$phi[,55]);
cov(flist$phi[,1],flist$phi[,55]);
cov(flist$phi[,2],flist$phi[,50]);

