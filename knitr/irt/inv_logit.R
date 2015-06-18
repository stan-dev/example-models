inv_logit <- function(u) {
  return(1 / (1 + exp(-u)));
}
