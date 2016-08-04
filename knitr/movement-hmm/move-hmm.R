library(moveHMM);

# downoad data with relevant headers for id, position, and movement
trackData <-
 read.table("http://www.esapubs.org/archive/ecol/E085/072/elk_data.txt",
            sep = "\t",
            header = TRUE)[1:735, c(1, 2, 3, 7)];

# id
lvls <- levels(trackData$Individual)[2:5];

sizes <- rep(NA, 4);
for (k in 1:4) sizes[k] <- sum(trackData$Individual == lvls[k]);
N1 <- sizes[1];  N2 <- sizes[2];  N3 <- sizes[3];  N4 <- sizes[4];
# identifies which elk
id <- c(rep(1, sizes[1]),
        rep(2, sizes[2]),
        rep(3, sizes[3]),
        rep(4, sizes[4]));

# movement (kilometers)
distance <- function(a1, a2) sqrt(sum((a1 - a2)^2));
dist <- rep(NA, 734);
for (n in 1:734) {
  dist[n] <- distance(c(trackData$Easting[n],
                        trackData$Northing[n]),
                      c(trackData$Easting[n + 1],
                        trackData$Northing[n + 1])) / 1000;
}

# bearings
bearing <- rep(0, 734);
delta_x <- rep(NA, 734);
delta_y <- rep(NA, 734);
for (n in 1:734) {
  delta_x[n] <- trackData$Easting[n + 1] - trackData$Easting[n];
  delta_y[n] <- trackData$Northing[n + 1] - trackData$Northing[n];
  if (delta_x[n] != 0 && delta_y[n] != 0) {
    bearing[n] <- atan2(delta_y[n], delta_x[n]);
  }
}

# change in bearings
# logic from turnAngle in moveHMM:
# https://github.com/cran/moveHMM/blob/master/R/turnAngle.R
turn <- rep(0, 735);
for (n in 2:734) {
  turn[n] <- bearing[n] - bearing[n - 1];
  while (turn[n] <= -pi) turn[n] <- turn[n] + 2 * pi;
  while (turn[n] > pi) turn[n] <- turn[n] -2 * pi;
}

# don't keep first step as no change in bearing (could use move)
# don't keep last step as no step from there
keep <- rep(TRUE, 735);
boundary <- 0;
for (k in 1:4) {
  keep[boundary + 1] <- FALSE;          # drop first point
  boundary <- boundary + sizes[k];
  keep[boundary] <- FALSE;          # drop last point
}

water <- trackData$dist_water..meters.[1:735] / 1000;

id <- id[keep];
dist <- dist[keep];
turn <- turn[keep];
x <- trackData$Easting[keep];
y <- trackData$Northing[keep];
water <- water[keep];

sizes <- sizes - 2;

df <- data.frame(id, dist, turn, x, y, water);

dist <- dist[1:sizes[1]];
turn <- turn[1:sizes[1]];
N <- length(dist);

library(rstan);
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
K <- 2;
m <- stan_model("move-hmm.stan");
fit <- sampling(m, data=c("N", "turn", "dist", "K"),
                chains = 1, iter = 200,
                control=list(stepsize=0.05, adapt_delta=0.9), refresh=2);
opt <- optimizing(m, data=c("N", "turn", "dist", "K"));
