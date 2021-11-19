data {
  int J;
  array[J] int n;
  vector[J] x;
  array[J] int y;
  real r;
  real R;
}
parameters {
  real<lower=0> sigma_theta;
  real<lower=0> sigma_distance;
  real<lower=0> sigma_y;
  real<lower=0> overshot;
  real<lower=0> distance_tolerance;
}
model {
  vector[J] p_angle;
  vector[J] p_distance;
  vector[J] p;
  p_angle = 2 * Phi(asin((R - r) ./ x) / sigma_theta) - 1;
  p_distance = Phi((distance_tolerance - overshot)
                   ./ ((x + overshot) * sigma_distance))
               - Phi((-overshot) ./ ((x + overshot) * sigma_distance));
  p = p_angle .* p_distance;
  to_vector(y) ./ to_vector(n) ~ normal(p,
                                        sqrt(p .* (1 - p) ./ to_vector(n)
                                             + sigma_y ^ 2));
  sigma_theta ~ normal(0, 1);
  sigma_distance ~ normal(0, 1);
  sigma_y ~ normal(0, 1);
  overshot ~ normal(1, 5);
  distance_tolerance ~ normal(3, 5);
}
generated quantities {
  real sigma_degrees;
  sigma_degrees = (180 / pi()) * sigma_theta;
}
