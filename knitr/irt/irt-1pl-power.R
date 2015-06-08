## ---- power-2pl ----
library(rstan);

## TEST 1a: 20 evenly spaced questions (-5 to 5)
b <- ((0:20) - 10) / 2;
J <- length(a);
a <- rep(1, J);



## TEST 1b: 20 evenly spaced questions (-5 to 5) discrim
a <- rep(4, J);

## TEST 1c: 20 evenly spaced (-5 to 5) bad discrim
a <_ rep(0.25, J);

## TEST 2a: 100 even spaced questions (-5 to 5)
b <- (((0:100) - 50) / 10;
J <- length(b);
a <- rep(1, J);




