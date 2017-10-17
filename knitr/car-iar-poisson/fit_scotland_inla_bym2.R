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

library(INLA)
# Set up the data
inla_data = list(
  y=data$y, E=data$E, x= 0.1 * data$x, region=c(1:data$N)
)

# Get the adjacency matrix
adj.matrix = sparseMatrix(i=nbs$node1,j=nbs$node2,x=1,symmetric=TRUE) 
adj.matrix = adj.matrix + Diagonal(nrow(adj.matrix))


# Set prior so that the mixing parameter 
# (rho in the Stan call, phi in INLA) has a beta(1/2,1/2) prior
# NB: INLA needs a prior on logit(phi)
# NB: INLA needs a prior for log(precision)
priors =  list(phi = list(prior="logitbeta",params=c(0.5,0.5)),
               prec = list(prior="logtgaussian",params=c(0,1/25))
               ) 

inla_formula <- y ~ 1+ x+ f(region, model = "bym2",graph=adj.matrix,hyper=priors,constr = TRUE)

inla_bym2 <- inla(inla_formula, family = "poisson", E=E, 
                  data = inla_data,
                  control.fixed = list(prec = 1/25,prec.intercept=1/25),
                  control.predictor = list(compute=TRUE), 
                  control.inla = list(strategy="laplace",fast=FALSE),
                  verbose=TRUE, debug=TRUE, silent=FALSE)

inla_bym2$summary.fixed
