library(rstan)   
options(mc.cores = parallel::detectCores())  

source("mungeCARdata4stan.R")  
source("scotland_data.R")
y = data$y;
x = 0.1 * data$x;
E = data$E;



nbs = mungeCARdata4stan(data$adj, data$num);
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

#Calculate the scaling of the model. Requires the INLA package.
#For convenience, this code isn't run.
if(FALSE) {
  library(INLA)
  #Build the adjacency matrix
  adj.matrix = sparseMatrix(i=nbs$node1,j=nbs$node2,x=1,symmetric=TRUE)
  #The ICAR precision matrix (note! This is singular)
  Q=  Diagonal(nbs$N, rowSums(adj.matrix)) - adj.matrix
  
  #Add a small jitter to the diagonal for numerical stability (optional but recommended)
  Q_pert = Q + Diagonal(nbs$N) * max(diag(Q)) * sqrt(.Machine$double.eps)


# Compute the diagonal elements of the covariance matrix subject to the 
# constraint that the entries of the ICAR sum to zero.
#See the function help for further details.
Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,nbs$N),e=0))

#Compute the geometric mean of the variances, which are on the diagonal of Q.inv
scaling_factor = exp(mean(log(diag(Q_inv))))

} else{
  scaling_factor= 0.4853175 # Not a magic number! Calculated as above
}

scot_stanfit = stan("bym2_predictor_plus_offset.stan",
         data=list(N,N_edges,node1,node2,y,x,E,scaling_factor),
         iter=10000);

print(scot_stanfit,
      pars=c("lp__", "beta0", "beta1", "rho", "sigma", "mu[5]"),
      probs=c(0.025, 0.5, 0.975));

