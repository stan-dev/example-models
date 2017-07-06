#-----------------------------------------------------------------------
# cgr
# Stan-based implementation of Cook, Gelman, Rubin 2006
# validate model by replication
#-----------------------------------------------------------------------
library(rstan)   
options(mc.cores = parallel::detectCores())  
rstan_options(auto_write=TRUE);



#-----------------------------------------------------------------------
# generate_data
# generate simulated data and params from model and fixed data
# return one draw from the simulation
# return null if cannot generate data from random inits
#
# params:
#  datagen_model compiled Stan model to generate data
#  datagen_data list of data inputs to model
#  num_chains int number of chains, default is 4
#  num_iters int total number of iterations, default is 2000
#            total draws will be num_chains * num_iters / 2
#            default number of draws is 4000
#  num_draw which draw to return, default is 4000
#
#  returns:
#   vector of values from num_draw-th draw
#-----------------------------------------------------------------------
generate_data <- function(datagen_model, datagen_data,
                          num_chains=4, num_iters=2000,
                          num_draw=4000) {
  total_draws = 0.5 * num_chains * num_iters;
  if (total_draws < num_draw) {
    print("cannot generate data, mismatch between specified chains, iterations, and draw");
    return(NULL);
  }  

  sim_data_fit = sampling(datagen_model, data=datagen_data,
                          chains=num_chains, iter=num_iters);
  
  sim_data_matrix = as.matrix(sim_data_fit);
  if (dim(sim_data_matrix)[1] != total_draws) {
    print("warning: simulation failed to generate data, frequent failures indicate problems with the data-generating Stan program");
    return(NULL);
  }

  result = sim_data_matrix[num_draw,];
  names(result) = colnames(sim_data_matrix);
  return(result);
}

#-----------------------------------------------------------------------
# get_stan_fit_quantiles
# fit model using simulated data and parameters (theta_0)
# compute the quantile of theta_0 w/r/t to samples from the posterior
# return vector of quantiles for all simulated parameters
#
# params:
#  draw vector of simulated values (created by generate_data)
#  param_names vector of names of parameters to check
#  target_model compiled Stan model to test
#  target_data list of simulated data inputs
#  num_chains int number of chains, default is 4
#  num_iters int total number of iterations, default is 2000
#
#  returns:
#   vector of percentile rank for each generated parameter
#-----------------------------------------------------------------------
get_stan_fit_quantiles <- function(draw, param_names, target_model, target_data,
                           num_chains=4, num_iters=2000) {
  target_fit = sampling(target_stan, data=target_data,
                        chains=num_chains, iter=num_iters);
  target_matrix = as.matrix(target_fit);

  rankings = rep(NA, length(param_names));
  mode(rankings) = "numeric";
  for (i in 1:length(param_names)) {
    theta0 = draw[param_names[i]];
    samples_theta0 = target_matrix[,param_names[i]];
    quantile = sum(samples_theta0 <= theta0)/length(samples_theta0);
    rankings[i] = quantile;
  }
  names(rankings) = param_names;
  return(rankings);
}

#-----------------------------------------------------------------------------
# main
#-----------------------------------------------------------------------------


DEBUG = TRUE;
num_replicates = 1000;

# compile models
datagen_model = stan_model("sim_bym_data.stan", model_name="sim_bym_model");
target_model = stan_model("bym_predictor_plus_offset.stan", model_name="bym_model");

# setup and initialize model-specific data
source("mungeCARdata4stan.R")  
source("scotland_data.R")
nbs = mungeCARdata4stan(data$adj, data$num);
E = data$E;
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;
fixed_data = list(N, N_edges, node1, node2, E)

num_chains = 4;
num_iters = 4000;
idx_draw = 8000;

ct_gen_fail = 0;  # track how often datagen fails

# generate first dataset
draw = NULL;
while (is.null(draw)) {
 draw = generate_data(datagen_model, fixed_data, num_chains, num_iters, idx_draw);
 if (is.null(draw)) ct_gen_fail = ct_gen_fail + 1;
}

# extract generated data for model from simulation from draw
data_names = names(draw);
x_idxs = grep("x\\[", data_names, value=FALSE);
y_idxs = grep("y\\[", data_names, value=FALSE);
x = as.vector(draw[x_idxs]);
y = as.vector(draw[y_idxs])
target_data = list(N, N_edges, node1, node2, E, y, x);

# create vector of names of params of interest
no_xs = setdiff(data_names, data_names[x_idxs]);
no_ys = setdiff(no_xs, data_names[y_idxs]);
param_names = setdiff(no_ys, c("lp__"));
std_param_names = grep(".*_std", param_names);
param_names = setdiff(param_names, param_names[std_param_names]);

# create matrix to record test replicates, one row per replicate
replicates_matrix = matrix(data=NA, nrow=0, ncol=length(param_names));
colnames(replicates_matrix) = param_names;

# do one test
quantiles = get_stan_fit_quantiles(draw, param_names, target_model, target_data, num_chains, num_iters);
replicates_matrix = rbind(replicates_matrix, quantiles);

# repeat procedure: data generation, model fit and quantile compuation 
for (i in 1:(num_replicates-1)) {
 if (DEBUG) {
   print(i);
   if (i %% 100 == 0) dput(replicates_matrix,file="tmp_snapshot.txt");
 }
 draw = NULL;
 while (is.null(draw)) {
   draw = generate_data(datagen_model, fixed_data, num_chains, num_iters, idx_draw);
   if (is.null(draw)) ct_gen_fail = ct_gen_fail + 1;
 }
 # update generated data for model
 x = as.vector(draw[x_idxs]);
 y = as.vector(draw[y_idxs])
 target_data = list(N, N_edges, node1, node2, E, y, x);

 quantiles = get_stan_fit_quantiles(draw, param_names, target_model, target_data, num_chains, num_iters);
 replicates_matrix = rbind(replicates_matrix, quantiles);
}
