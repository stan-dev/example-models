## template of an R script to fit BYM2 model

library(rstan);
options(mc.cores = 2);

library(maptools);   
library(spdep);
library(rgdal)

# using nb_object helper functions
source("nb_data_funs.R");

# load data - need outcome, exposure
#y = mydata...
#E = mydata...

# if co-variates, set up K and x
#K = 0;
#x = matrix(nrow=0,ncol=0);

# load the fully connected nb object
# or load shapefiles and do `poly2nb`
#nb = ...

## bym2 model requires symmetric neighborhood structure
validateNb(nb);

## check that nb object is fully connected.
!(isDisconnected(nb))

# calculate scaling factor
scaling_factor = scale_nb_components(nb)[1];

# entire graph represented as node pairs
nbs=nb2graph(nb);
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;

# fit model
bym2_stanfit = stan("bym2.stan",
  data=list(N,N_edges,node1,node2,y,E,scaling_factor),
  iter=5000, warmup=4000, chains= 4, control = list(adapt_delta = 0.95), save_warmup=FALSE);
