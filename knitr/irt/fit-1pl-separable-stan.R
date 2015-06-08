y_sep <- matrix(NA, I + 1, J);
y_sep[1:I, 1:J] <- y
I_sep <- I + 1
y_sep[I_sep,] <- rep(0, J);

fit <- stan("irt_1pl_mle.stan", data = list(I=I_sep, J=J, y=y_sep));
