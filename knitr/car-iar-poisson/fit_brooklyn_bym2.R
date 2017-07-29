library(rstan);
options(mc.cores = parallel::detectCores());
library(maptools);   
library(spdep);

source("nyc_ped_subset.data.R");
y = events_all_2001[all_tractIDs %in% bklyn_tractIDs];
x = pop_adj_2001[all_tractIDs %in% bklyn_tractIDs];

source("nbdata4stan.R");
nyc_all_tracts.shp<-readShapePoly("nycTracts10/nycTracts10");
bklyn_tracts <- nyc_all_tracts.shp$GEOID10 %in% bklyn_tractIDs;
bklyn_tracts.shp <- nyc_all_tracts.shp[bklyn_tracts,];
bklyn_tracts.shp <- bklyn_tracts.shp[order(bklyn_tracts.shp$GEOID10),];
nb_bk = poly2nb(bklyn_tracts.shp);
nbs=nbdata4stan(nb_bk);
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

DO_CALC = TRUE;
#Calculate the scaling of the model. Requires the INLA package.
if(DO_CALC) {
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
  scaling_factor= 0.5;
}

bym2_bk_stanfit = stan("bym2_predictor_only.stan", data=list(N,N_edges,node1,node2,y,x,scaling_factor), iter=6000);
print(bym2_bk_stanfit, digits=3, pars=c("lp__", "beta0", "beta1", "rho", "sigma", "mu[1]", "mu[2]", "mu[3]", "mu[701]", "mu[702]", "mu[703]"), probs=c(0.025, 0.5, 0.975));

