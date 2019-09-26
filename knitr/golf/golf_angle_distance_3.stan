data {
  int J;
  int n[J];
  vector[J] x;
  int y[J];
  real r;
  real R;
  real overshot;
  real distance_tolerance;
}
transformed data {
  vector[J] threshold_angle = asin((R-r) ./ x);
  vector[J] raw_proportions = to_vector(y) ./ to_vector(n);
}
parameters {
  real<lower=0> sigma_angle;
  real<lower=0> sigma_distance;
  real<lower=0> sigma_y;
}
model {
  vector[J] p_angle = 2*Phi(threshold_angle / sigma_angle) - 1;
  vector[J] p_distance = Phi((distance_tolerance - overshot) ./ ((x + overshot)*sigma_distance)) -
               Phi((- overshot) ./ ((x + overshot)*sigma_distance));
  vector[J] p = p_angle .* p_distance;
  raw_proportions ~ normal(p, sqrt(p .* (1-p) ./ to_vector(n) + sigma_y^2));
  [sigma_angle, sigma_distance, sigma_y] ~ normal(0, 1);
}
generated quantities {
  real sigma_degrees = sigma_angle * 180 / pi();
}
