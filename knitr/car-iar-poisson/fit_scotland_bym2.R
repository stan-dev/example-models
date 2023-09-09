library(devtools)
if(!require(cmdstanr)){
  devtools::install_github("stan-dev/cmdstanr", dependencies=c("Depends", "Imports"))
}
if(!require(INLA)){
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}
library(cmdstanr)   
library(INLA)

source("mungeCARdata4stan.R")  
source("scotland_data.R")
y = data$y;
x = 0.1 * data$x;
E = data$E;
K = 1;

nbs = mungeCARdata4stan(data$adj, data$num);
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

#Build the adjacency matrix using INLA library functions
adj.matrix = sparseMatrix(i=nbs$node1,j=nbs$node2,x=1,symmetric=TRUE)
#The ICAR precision matrix (note! This is singular)
Q=  Diagonal(nbs$N, rowSums(adj.matrix)) - adj.matrix
#Add a small jitter to the diagonal for numerical stability (optional but recommended)
Q_pert = Q + Diagonal(nbs$N) * max(diag(Q)) * sqrt(.Machine$double.eps)

# Compute the diagonal elements of the covariance matrix subject to the 
# constraint that the entries of the ICAR sum to zero.
#See the inla.qinv function help for further details.
Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,nbs$N),e=0))

#Compute the geometric mean of the variances, which are on the diagonal of Q.inv
scaling_factor = exp(mean(log(diag(Q_inv))))

data = list(N=N,
            N_edges=N_edges,
            node1=node1,
            node2=node2,
            y=y,
            x=x,
            E=E,
            scaling_factor=scaling_factor);

bym2_model = cmdstan_model("bym2_predictor_plus_offset.stan");

bym2_scot_stanfit = bym2_model$sample(
                                   data=data,
                                   parallel_chains=4,
                                   refresh=0);

bym2_scot_stanfit$summary(variables = c("beta0", "beta1",
                                       "sigma", "rho",
                                       "mu[5]","phi[5]","theta[5]"))

