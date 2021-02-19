library(SpatialEpi)
data("scotland")

source("scotland_nbs.data.R")
source("scotland_3_comp_nbs.data.R")
source("scotland_islands_nbs.data.R")


# plots
plot(scot_nb, scotland$spatial.polygon)
plot(scot_3_comp_nb, scotland$spatial.polygon)
plot(scot_islands_nb, scotland$spatial.polygon)

source("bym2_helpers.R")

# connected graph
scot_nodes = length(scot_nb)
scot_edge_array = nb_to_edge_array(scot_nb)
scot_edges = dim(scot_edge_array)[2]
scot_scaling_factor = scaling_factor(scot_edge_array)

# disconnected graph - 2 components, 1 singleton (island)

scot_3_comp_nodes = length(scot_3_comp_nb)
scot_3_comp_edge_array = nb_to_edge_array(scot_3_comp_nb)
scot_3_comp_etc = index_components(scot_3_comp_nb)

# disconnected graph - 1 components, 3 singletons (island)

scot_islands_nodes = length(scot_islands_nb)
scot_islands_edge_array = nb_to_edge_array(scot_islands_nb)
scot_islands_etc = index_components(scot_islands_nb)
