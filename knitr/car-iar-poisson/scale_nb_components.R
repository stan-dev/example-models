library(INLA);
library(spdep);

# indexByComponent
#
# input: vector of component ids
# returns: vector of per-component consecutive node ids
#
renumber = function(x) {
  y = x;
  comps = as.matrix(table(x));
  num_comps = nrow(comps);
  for (i in 1:nrow(comps)) {
    idx = 1;
    rel_idx = 1;
    while (idx <= length(x)) {
      if (x[idx] == i) {
        y[idx] = rel_idx;
        rel_idx = rel_idx + 1;
      }
      idx = idx + 1;
    }
  }
  return(y);
}

# for a given subcomponent, return graph as lists of node1, node2 pairs
#
# inputs:
# x: nb object
# c_id: subcomponent id
# comp_ids: vector of subcomponent ids
# offsets: vector of subcomponent node numberings
# returns: list of node1, node2 ids
#
nb2subgraph = function(x, c_id, comp_ids, offsets) {
  N = length(x);
  n_links = 0;
  for (i in 1:N) {
    if (comp_ids[i] == c_id) {
      if (x[[i]][1] != 0) {
        n_links = n_links + length(x[i]);
      }
    }
  }
  N_edges = n_links / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  idx = 0;
  for (i in 1:N) {
    if (comp_ids[i] == c_id) {
      if (x[[i]][1] != 0) {
        for (j in 1:length(x[[i]])) {
          n2 = unlist(x[[i]][j]);    
          if (i < n2) {
            idx = idx + 1;
            node1[idx] = offsets[i];
            node2[idx] = offsets[n2];
          }
        }
      }
    }
  }
  return (list("node1"=node1,"node2"=node2));
}

# scale_nb_components
#
# input: nb_object
# returns: vector of length(nb_object) where each element is scaling factor for component 
# scaling factor for singletons is 0
#
scale_nb_components = function(x) {
  N = length(x);
  scales = rep(0, N);
  comp_ids = n.comp.nb(x)[[2]];
  offsets = indexByComponent(comp_ids);
  comps = as.matrix(table(comp_ids));
  num_comps = nrow(comps);
  for (i in 1:num_comps) {
    N_subregions = comps[i,1];
    if (N_subregions > 1) {
      # get adj matrix for this component
      drops = comp_ids != i;
      nb_tmp = droplinks(x, drops);      
      nb_graph = nb2subgraph(nb_tmp, i, comp_ids, offsets);
      adj.matrix = sparseMatrix( i=nb_graph$node1, j=nb_graph$node2, x=1, dims=c(N_subregions,N_subregions), symmetric=TRUE);
      # compute ICAR precision matrix
      Q =  Diagonal(N_subregions, rowSums(adj.matrix)) - adj.matrix;
      # Add a small jitter to the diagonal for numerical stability (optional but recommended)
      Q_pert = Q + Diagonal(N_subregions) * max(diag(Q)) * sqrt(.Machine$double.eps)
      # Compute the diagonal elements of the covariance matrix subject to the 
      # constraint that the entries of the ICAR sum to zero.
      Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,N_subregions),e=0))
      # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
      scaling_factor = exp(mean(log(diag(Q_inv))))

      scales[comp_ids == i] = scaling_factor;
    } else {
      scales[comp_ids == i] = -1;
    }
  }
  return(scales);
}

# nb2graph
#
# input: nb_object
# returns: dataframe containing num nodes, num edges,
#          and a list of graph edges from node1 to node2.
#
# note: Stan requires consecutive node numbers because there's a 1:1
# correspondance between nodes and params in spatial model;
# must skip island nodes and adjust subsequenct indices accordingly.
#
nb2graph = function(x) {
  N = length(x);
  n_links = 0;
  for (i in 1:N) {
    if (x[[i]][1] != 0) {
      n_links = n_links + length(x[[i]]);
    }
  }
  N_edges = n_links / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  idx = 0;
  adjust = 0;
  for (i in 1:N) {
    if (x[[i]][1] == 0) {
      adjust = adjust + 1;
    } else {
      for (j in 1:length(x[[i]])) {
        n2 = unlist(x[[i]][j]);
        if (i < n2) {
          idx = idx + 1;
          node1[idx] = i - adjust;
          node2[idx] = n2 - adjust;
        }
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}

