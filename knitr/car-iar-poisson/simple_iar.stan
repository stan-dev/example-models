data {
  int<lower=0> N;                 // num regions
  int<lower=0> N_edges;           // num edges, (undirected)
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // node1[i] adjacent to node2[i]
}
parameters {
  vector[N-1] phi_raw_std;
}
transformed parameters {
  vector[N] phi;
  phi[1:(N - 1)] = phi_raw_std;
  phi[N] = -sum(phi_raw_std);
}
model {
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
}
