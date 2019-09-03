library(INLA);
library(spdep);


# validateNb
# check that nbObject is symmetric, has connected components
#
validateNb = function(x) {
  if (is.symmetric.nb(x) && n.comp.nb(x)[[1]] < length(x)) return(TRUE);
  return(FALSE);
}


# isDisconnected
# check that nbObject has more than 1 component
#
isDisconnected = function(x) {
  return(n.comp.nb(x)[[1]] > 1);
}


# indexByComponent
#
# input: vector of component ids
# returns: vector of per-component consecutive node ids
#
indexByComponent = function(x) {
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


# nb2graph
#
# input: nb_object
# returns: dataframe containing num nodes, num edges,
#          and a list of graph edges from node1 to node2.
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
  for (i in 1:N) {
    if (x[[i]][1] > 0) {
      for (j in 1:length(x[[i]])) {
        n2 = unlist(x[[i]][j]);
        if (i < n2) {
          idx = idx + 1;
          node1[idx] = i;
          node2[idx] = n2;
        }
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}


# nb2subgraph
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


# orderByComp
# given nbObject, reorder nb object so that all components are contiguous
# singletons moved to end of list
# returns list containing:
#   - new nbObject
#   - vector which maps old index to new index
#
orderByComponent = function(x) {
  if (!validateNb(x)) return(list());
  N = length(x);
  if (!isDisconnected(x)) return(list(x,seq(1:N)));

  rMap = rep(integer(0),N);
  comp_ids = n.comp.nb(x)[[2]];
  comps = as.matrix(table(comp_ids));
  num_comps = nrow(comps);
  comp_sizes = as.vector(comps[,1]);
  idx = 1;
  for (i in 1:nrow(comps)) {
    if (comp_sizes[i] > 1) {
      positions = which(comp_ids == i);
      for (j in 1:length(positions)) {
        rMap[idx] = as.integer(positions[j])
        idx = idx + 1;
      }
    }
  }
  for (i in 1:nrow(comps)) {
    if (comp_sizes[i] == 1) {
      positions = which(comp_ids == i);
      for (j in 1:length(positions)) {
        rMap[idx] = as.integer(positions[j])
        idx = idx + 1;
      }
    }
  }
  new_ids = vector("character", length=N);
  for (i in 1:N) {
    idx_old = rMap[i];
    new_ids[i] = attributes(x)$region.id[idx_old];
  }

  # generate new nb list
  new_nb = structure(vector("list", length=N),class="nb");
  attr(new_nb, "region.id") = new_ids;  
  attr(new_nb, "type") = attributes(x)$type;
  attr(new_nb, "sym") = attributes(x)$sym;
  attr(new_nb, "region.id") = new_ids;  
  for (i in 1:N) {
    idx_old = rMap[i];
    old_nbs = x[[idx_old]];
    num_nbs = length(old_nbs);
    new_nb[[i]] = vector("integer", length=num_nbs);
    for (j in 1:num_nbs) {
      old_id = old_nbs[j];
      if (old_id == 0) {
        new_nb[[i]][j] = as.integer(0);
      } else {
        new_id = which(rMap == old_id);
        new_nb[[i]][j] = as.integer(new_id);
      }
    }
  }
  return(list(new_nb,rMap));
}  


# reorderVector
#
# input: data vector, offsets vector
# returns: vector of same length as input data vector
#          reordered according to offsets
#
reorderVector = function(x, rMap) {
  if (!is.vector(x)) return(NULL);
  N = length(x);              
  result = vector("numeric", length=N);
  for (i in 1:N) {
    result[i]= x[rMap[i]];
  }
  return(result);
}


# reorderMatrix
#
# input: data matrix, offsets vector
# returns: matrix of same shape as input data matrix,
#          rows reordered according to offsets
#
reorderMatrix = function(x, rMap) {
  if (!is.matrix(x)) return(NULL);
  N = nrow(x);
  result = matrix(nrow=N, ncol=ncol(x));
  for (i in 1:N) {
    result[i,]= x[rMap[i],];
  }
  return(result);
}


# scale_nb_components
#
# input: nb_object
# returns: vector of per-component scaling factor (for BYM2 model)
# scaling factor for singletons is 0
#
scale_nb_components = function(x) {
  N = length(x);
  comp_ids = n.comp.nb(x)[[2]];
  offsets = indexByComponent(comp_ids);

  comps = as.matrix(table(comp_ids));
  num_comps = nrow(comps);
  scales = vector("numeric", length=num_comps);
  for (i in 1:num_comps) {
    N_subregions = comps[i,1];
    scales[i] = 0.0;
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
      scales[i] = scaling_factor;
    }
  }
  return(scales);
}

library(INLA)
get_scaling_factor = function(nbs) {
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
    return((mean(log(diag(Q_inv)))))
}
