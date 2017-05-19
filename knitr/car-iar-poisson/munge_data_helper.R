get_nb_off_diags = function(x) {
     n_tracts = length(x);
     n_links = sum(card(x));
     idx_i = vector(mode="numeric",length=n_links);
     idx_j = vector(mode="numeric",length=n_links);
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
     return(cbind(idx_i,idx_j));
}


  
  
  
