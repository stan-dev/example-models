functions {
  /**
   * Return the log probability density of the specified vector of
   * coefficients under the ICAR model with unit variance, where
   * adjacency is determined by the adjacency array and the spatial
   * structure is a disconnected graph which has at least one
   * connected component.  Each connected component has a
   * soft sum-to-zero constraint.   Singleton nodes have
   * distribution normal(0, 1/sqrt(K))
   *
   * The spatial structure is described by a 2-D adjacency array
   * over the all edges in the areal map and a arrays of the
   * indices of per-component nodes and edges which are used as
   * masks into phi and the adjacency matrix.   Because the Stan
   * language lacks ragged arrays, these are all square matrices,
   * padded out with zeros; additional vectors record the number
   * of nodes and edges in each component.
   *
   * @param phi vector of varying effects
   * @param adjacency parallel arrays of indexes of adjacent elements of phi
   * @param node_cts array of sizes of per-component nodes
   * @param edge_cts array of sizes of per-component edges
   * @param node_idxs array of arrays of per_component node indexes.
   * @param edge_idxs array of arrays of per_component edge indexes.
   *
   * @return ICAR log probability density
   *
   * @reject if the the adjacency matrix does not have two rows
   * @reject if size mismatch between indexing arrays
   * @reject if size mismatch between phi and node indexes columns.
   */
  real standard_icar_disconnected_lpdf(vector phi, array[,] int adjacency,
                                       array[] int node_cts,
                                       array[] int edge_cts,
                                       array[,] int node_idxs,
                                       array[,] int edge_idxs) {
    int num_nodes = size(phi);
    int num_edges = dims(adjacency)[2];
    int num_comps = size(edge_cts);
    if (size(adjacency) != 2) {
      reject("require 2 rows for adjacency array;", " found rows = ",
             size(adjacency));
    }
    if (!(num_nodes == dims(node_idxs)[2] && size(node_cts) == size(edge_cts)
          && size(node_cts) == size(node_idxs)
          && size(edge_cts) == size(edge_idxs))) {
      reject("arguments have size mismatch, expecting ", num_comps,
             " rows for node_cts edge_cts, node_idxs, and edge_idxs,",
             num_nodes, " elements in phi and columns of node_idxs, and ",
             num_edges, " columns of edge_idxs.");
    }
    
    real total = 0;
    for (n in 1 : num_comps) {
      if (node_cts[n] > 1) {
        total += -0.5
                 * dot_self(phi[adjacency[1, edge_idxs[n, 1 : edge_cts[n]]]]
                            - phi[adjacency[2, edge_idxs[n, 1 : edge_cts[n]]]])
                 + normal_lpdf(sum(phi[node_idxs[n, 1 : node_cts[n]]]) | 0, 0.001
                                                                    * node_cts[n]);
      } else {
        total += normal_lpdf(phi[node_idxs[n, 1]] | 0, 1);
      }
    }
    return total;
  }
}
data {
  // spatial structure
  int<lower=0> I; // number of nodes
  int<lower=0> J; // number of edges
  array[2, J] int<lower=1, upper=I> edges; // node[1, j] adjacent to node[2, j]
  
  int<lower=0, upper=I> K; // number of components in spatial graph
  array[K] int<lower=0, upper=I> K_node_cts; // per-component nodes
  array[K] int<lower=0, upper=J> K_edge_cts; // per-component edges
  array[K, I] int<lower=0, upper=I> K_node_idxs; // rows contain per-component node indexes
  array[K, J] int<lower=0, upper=J> K_edge_idxs; // rows contain per-component edge indexes
  
  vector[K] tau; // scaling factor
  
  array[I] int<lower=0> y; // count outcomes
  vector<lower=0>[I] E; // exposure
  vector[I] x; // predictor
}
transformed data {
  vector[I] log_E = log(E);
}
parameters {
  real alpha; // intercept
  real beta; // covariates
  
  // spatial effects
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance
  real<lower=0> sigma; // scale of spatial effects
  vector[I] theta; // standardized heterogeneous spatial effects
  vector[I] phi; // standardized spatially smoothed spatial effects
}
transformed parameters {
  vector[I] gamma;
  for (k in 1 : K) {
    gamma[K_node_idxs[k, 1 : K_node_cts[k]]] = (sqrt(1 - rho)
                                                * theta[K_node_idxs[k, 1 : K_node_cts[k]]]
                                                + sqrt(rho / tau[k])
                                                  * phi[K_node_idxs[k, 1 : K_node_cts[k]]])
                                               * sigma;
  }
}
model {
  y ~ poisson_log(log_E + alpha + x * beta + gamma * sigma); // co-variates
  
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  
  // spatial hyperpriors and priors
  sigma ~ normal(0, 1);
  rho ~ beta(0.5, 0.5);
  theta ~ normal(0, 1);
  phi ~ standard_icar_disconnected(edges, K_node_cts, K_edge_cts,
                                   K_node_idxs, K_edge_idxs);
}
generated quantities {
  // posterior predictive checks
  vector[I] eta = log_E + alpha + x * beta + gamma * sigma;
  vector[I] y_prime = exp(eta);
  //   int y_rep[I,10];
  //   for (j in 1:10) {
  //     if (max(eta) > 20) {
  //       // avoid overflow in poisson_log_rng
  //       print("max eta too big: ", max(eta));
  //       for (i in 1:I)
  //   	y_rep[i,j] = -1;
  //     } else {
  //       for (i in 1:I)
  //         y_rep[i,j] = poisson_log_rng(eta[i]);
  //     }
  //   }
  real logit_rho = log(rho / (1.0 - rho));
}
