get_nb_off_diags = function(x) {
     n_tracts = length(x);
     n_links = sum(card(x));
     idx_i = vector(mode="numeric", length=n_links);
     idx_j = vector(mode="numeric", length=n_links);
     idx = 1;
     for (i in 1:n_tracts) {
        if (x[[i]][1] != 0) {
           for (j in 1:length(x[[i]])) {
              idx_i[idx] = i;
              idx_j[idx] = x[[i]][j];
              idx = idx + 1;
           }
        }
     }
     t1 =  cbind(idx_i, idx_j);
     t2 = aggregate(count~idx_i+idx_j, transform(t1, idx_i=pmin(idx_i, idx_j), idx_j=pmax(idx_i, idx_j), count=1), sum)
     t3 = data.frame(t2$idx_i, t2$idx_j);
     return (as.matrix(t3));
}
