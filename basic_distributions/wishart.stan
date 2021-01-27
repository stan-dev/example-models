transformed data {
  cov_matrix[2] S = [[2, 0],
                     [0, 0.5]];
} 
parameters {
  cov_matrix[2] W; 
} 
model {
  W ~ wishart(4, S); 
} 
