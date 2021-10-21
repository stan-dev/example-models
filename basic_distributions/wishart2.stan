transformed data {
  cov_matrix[4] S;
  
  S[1, 1] = 2.9983662;
  S[2, 1] = 0.2898776;
  S[3, 1] = -2.650523;
  S[4, 1] = 0.1055911;
  S[2, 2] = 11.4803610;
  S[3, 2] = 7.157993;
  S[4, 2] = -3.1129955;
  S[3, 3] = 11.676181;
  S[4, 3] = -3.5866852;
  S[4, 4] = 1.4482736;
  
  S = symmetrize_from_lower_tri(S);
}
parameters {
  cov_matrix[4] W;
}
model {
  W ~ wishart(10, S);
}
