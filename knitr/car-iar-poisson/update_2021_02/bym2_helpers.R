library(Matrix)
library(INLA);
library(igraph)
library(sf)
library(spdep)

## Helper functions which take R package spdep nb objects
## to input data needed by Stan implementation of BYM2 model.
##
## Functions compute:
##
##  - spatial structure as edge array
##  - scaling factor for spatial, heterogeneous effects in BYM2 model.
##  - for graphs which disconnected components, need to identify
##    per-component edges, scaling factors

# convert nb object to array of column-wise edge pairs
# creates a 2 X N array where N is the number of edges
# and entry [n,1] and [n,2] contain the index of the nodes
# connected by edge n.
nb_to_edge_array <- function(nb_obj) {
    adj_matrix = nb2mat(nb_obj,style="B", zero.policy = TRUE)
    t(as_edgelist(graph_from_adjacency_matrix(adj_matrix, mode="undirected")))
}


# compute geometric mean of a vector
geometric_mean <- function(x) exp(mean(log(x))) 


# compute scaling factor for a fully connected areal map
# accounts for differences in spatial connectivity
scaling_factor <- function(edge_array) {
    adj_matrix = sparseMatrix(i=edge_array[1, ],j=edge_array[2, ],x=1,symmetric=TRUE)

    N = dim(adj_matrix)[1]

    # Create ICAR precision matrix  (diag - adjacency): this is singular
    # function Diagonal creates a square matrix with given diagonal
    Q =  Diagonal(N, rowSums(adj_matrix)) - adj_matrix

    # Add a small jitter to the diagonal for numerical stability (optional but recommended)
    Q_pert = Q + Diagonal(N) * max(diag(Q)) * sqrt(.Machine$double.eps)

    # Function inla.qinv provides efficient way to calculate the elements of the
    # the inverse corresponding to the non-zero elements of Q
    Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,N),e=0))

    # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
    scaling_factor <- geometric_mean(Matrix::diag(Q_inv)) 
    return(scaling_factor) 
}


# for disconnected graph, return list of comp sizes, and 2-D arrays
# which index nodes and edges per component, respectively.
index_components <- function(nb_obj) {
    num_nodes = length(nb_obj)

    comps = n.comp.nb(nb_obj)
    comp_idxs = comps$comp.id
    num_comps = comps$nc
    all_edges = nb_to_edge_array(nb_obj)
    num_edges = dim(all_edges)[2]

    comp_node_idxs = matrix(data=c(0), nrow=num_comps, ncol=num_nodes)
    for (k in 1:num_comps) {
        comp_node_idx = which(comp_idxs == k)
        while (length(comp_node_idx) < num_nodes) {
            comp_node_idx[(length(comp_node_idx)+1)] = 0
        }
        comp_node_idxs[k,] = comp_node_idx
    }

    edges_per_comp = vector(mode="integer", length=num_comps)
    comp_edge_idxs = matrix(data=c(0), nrow=num_comps, ncol=num_edges)
    for (k in 1:num_comps) {
        comp_idx = which(comp_idxs == k)
        comp_edge_idx = which(all_edges[1,] %in% comp_idx)
        edges_per_comp[k] = length(comp_edge_idx)
        while (length(comp_edge_idx) < num_edges) {
            comp_edge_idx[(length(comp_edge_idx)+1)] = 0
        }
        comp_edge_idxs[k,] = comp_edge_idx
    }

    scaling_factors = vector(mode="numeric", length=num_comps)
    for (k in 1:num_comps) {
        if (edges_per_comp[k] > 0) {
            sub_nb = subset(nb_obj, comp_idxs == k)
            edges = nb_to_edge_array(sub_nb)
            scaling_factors[k] = scaling_factor(edges)
        } else {
            scaling_factors[k] = 1.0
        }
    }   

    return(list("K"=num_comps,
                "K_num_nodes"=as.vector(table(comp_idxs)),
                "K_num_edges"=edges_per_comp,
                "K_node_idxs"=comp_node_idxs,
                "K_edge_idxs"=comp_edge_idxs,
                "K_scaling_factors"=scaling_factors))
}

