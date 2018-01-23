library(rstan)   
options(mc.cores = parallel::detectCores())  

source("mungeCARdata4stan.R")  
source("scotland_data.R")

nbs = mungeCARdata4stan(data$adj, data$num);
N = data$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

iar_stan = stan_model("simple_iar.stan");
fit_stan = sampling(iar_stan, data=list(N,N_edges,node1,node2), control=list(adapt_delta = 0.97, stepsize = 0.1), chains=2, warmup=9000, iter=10000, save_warmup=FALSE);

ofile = "simple_iar_stan_010K_iters.txt";
capture.output(print(fit_stan, digits=3, probs=c(0.025, 0.975)),file=ofile);


