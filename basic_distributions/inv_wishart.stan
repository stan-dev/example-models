transformed data {
  cov_matrix[3] S = [[2, 0, 0], [0, 1, 0], [0, 0, 0.5]];
}
parameters {
  cov_matrix[3] W;
}
model {
  W ~ inv_wishart(5, S);
}
