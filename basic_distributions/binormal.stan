transformed data {
    matrix[2,2] Sigma = [[1, 0.1],
                         [0.1, 1]];
    vector[2] mu = [0, 0]';
}
parameters {
    vector[2] y;
}
model {
      y ~ multi_normal(mu,Sigma);
}
