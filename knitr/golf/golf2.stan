data {
  int J;
  array[J] int n;
  vector[J] x;
  array[J] int y;
  real r;
  real R;
  real overshot;
  real distance_tolerance;
}
parameters {
  real<lower=0> sigma_theta;
  real<lower=0> sigma_distance;
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
  y ~ binomial(n, p);
}
generated quantities {
  real sigma_degrees;
  sigma_degrees = (180 / pi()) * sigma_theta;
}
