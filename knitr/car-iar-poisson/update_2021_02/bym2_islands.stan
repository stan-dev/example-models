functions {
  /**
   * Return the log probability density of the specified vector of
   * coefficients under the ICAR model with unit variance, where
   * adjacency is determined by the adjacency array and the spatial
   * structure is a disconnected graph which has at least one
   * connected component.  The spatial structure is described by
   * a 2-D adjacency array over the all edges in the areal map
   * and a series of arrays of per-component indexes.
   * The node arrays provide per-component masks into phi, the
   * per-node ICAR component, and corresponding masks into the
   * edge array.  Because the Stan language lacks ragged arrays,
   * these are all zero-padded square matrices, and additional vectors
   * record the number of nodes and edges in each component.
   *
   *
   * Each connected component has a soft sum-to-zero constraint.
   * Singleton components don't contribute to the ICAR model.
   *
   * @param phi vector of varying effects
   * @param adjacency parallel arrays of indexes of adjacent elements of phi
   * @param num_nodes array of sizes of per-component nodes
   * @param num_edges array of sizes of per-component edges
   * @param node_idxs array of arrays of per_component node indexes.
   * @param edge_idxs array of arrays of per_component edge indexes.
   *
   * @return ICAR log probability density
   *
   * @reject if the the adjacency matrix does not have two rows
   * @reject if size mismatch between indexing arrays
   * @reject if size mismatch between phi and dimension 2 of comp_members
   */
  real standard_icar_disconnected_lpdf(vector phi,
				       int[ , ] adjacency,
				       int[ ] num_nodes,
				       int[ ] num_edges,
				       int[ , ] node_idxs,
				       int[ , ] edge_idxs) {
    if (size(adjacency) != 2)
      reject("require 2 rows for adjacency array;",
             " found rows = ", size(adjacency));
      
    if (!(size(num_nodes) == size(num_edges)
	  && size(num_nodes) == size(node_idxs)
	  && size(num_edges) == size(edge_idxs)))
      reject("bad graph indexes, expecting ",
	     size(num_nodes),
	     " rows for node and edge index matrices;",
             " inputs have ",
	     size(num_nodes),
	     " and ",
	     size(num_edges),
	     " respectively.");

    real total = 0;
    for (n in 1:size(num_nodes)) {
      if (num_nodes[n] > 1)
	total += -0.5 * dot_self(phi[adjacency[1, edge_idxs[n, 1:num_edges[n]]]] -
				 phi[adjacency[2, edge_idxs[n, 1:num_edges[n]]]])
	  + normal_lpdf(sum(phi[node_idxs[n, 1:num_nodes[n]]]) | 0, 0.001 * num_nodes[n]);
      else
	total += normal_lpdf(phi[node_idxs[n, 1]] | 0, inv_sqrt(1.0 * size(num_nodes)));
    }
    return total;
  }
}
data {
  // spatial structure
  int<lower = 0> I;  // number of nodes
  int<lower = 0> J;  // number of edges
  int<lower = 1, upper = I> edges[2, J];  // node[1, j] adjacent to node[2, j]

  int<lower=0, upper=I> K;  // number of components in spatial graph
  int<lower=0, upper=I> K_num_nodes[K];   // per-component nodes
  int<lower=0, upper=J> K_num_edges[K];   // per-component edges
  int<lower=0, upper=I> K_node_idxs[K, I];  // rows contain per-component node indexes
  int<lower=0, upper=J> K_edge_idxs[K, J];  // rows contain per-component edge indexes

  vector[K] tau; // scaling factor

  int<lower=0> y[I];              // count outcomes
  vector<lower=0>[I] E;           // exposure
  vector[I] x;                 // predictor
}
transformed data {
  vector[I] log_E = log(E);
}
parameters {
  real alpha;            // intercept
  real beta;       // covariates

  // spatial effects
  real<lower = 0> sigma;  // scale of spatial effects
  real<lower = 0, upper = 1> rho;  // proportion of spatial effect that's spatially smoothed
  vector[I] theta;  // standardized heterogeneous spatial effects
  vector[I] phi;  // standardized spatially smoothed spatial effects
}
transformed parameters {
  vector[I] gamma;
  // each component has its own spatial smoothing
  for (k in 1:K) {
    if (K_num_nodes[k] == 1) {
      gamma[K_node_idxs[k,1]] =
	theta[K_node_idxs[k,1]] + normal_lpdf(phi[K_node_idxs[k,1]] | 0, inv_sqrt(K));
    } else {
      gamma[K_node_idxs[k, 1:K_num_nodes[k]]] = 
	    (sqrt(1 - rho) * theta[K_node_idxs[k, 1:K_num_nodes[k]]]
	     + (sqrt(rho) * sqrt(1 / tau[k])
		* phi[K_node_idxs[k, 1:K_num_nodes[k]]]) * sigma);
    }
  }
}
model {
  y ~ poisson_log(log_E + alpha + x * beta + gamma * sigma);  // co-variates

  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);

  // spatial hyperpriors and priors
  sigma ~ normal(0, 1);
  rho ~ beta(0.5, 0.5);
  theta ~ normal(0, 1);
  phi ~ standard_icar_disconnected(edges, K_num_nodes, K_num_edges, K_node_idxs, K_edge_idxs);
}
