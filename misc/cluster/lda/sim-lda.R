library("gtools");  # for rdirichlet

V <- 5; # words: river, stream, bank, money, loan

K <- 2; # topics: RIVER, BANK

phi <- array(NA,c(2,5));
phi[1,] = c(0.330, 0.330, 0.330, 0.005, 0.005);
phi[2,] = c(0.005, 0.005, 0.330, 0.330, 0.330);

set.seed(123)

M <- 25;  # docs
avg_doc_length <- 10;
doc_length <- rpois(M,avg_doc_length);
N <- sum(doc_length);

alpha <- rep(1/K,K);
beta <- rep(1/V,V);

theta <- rdirichlet(M,alpha);

corpus <- matrix(0, nrow = M, ncol = V);
for (m in 1:M) {
  for (i in 1:doc_length[m]) {
    topic_id <- which(rmultinom(1,1,theta[m,]) == 1);
    word_id <- which(rmultinom(1,1,phi[topic_id,]) == 1);
    corpus[m, word_id] = corpus[m, word_id] + 1;
  }
}
stopifnot(all(rowSums(corpus) == doc_length));

dump(c("K","V","M","corpus","alpha","beta"),"lda.data.R");