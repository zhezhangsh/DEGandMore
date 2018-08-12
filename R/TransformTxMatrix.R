# Transform matrix of transcriptome data using different methods
TransformTxMatrix <- function(d, g, method='rank2normal', rank=NA) {
  
  if (tolower(method[1])=='rank2normal' & length(rank)==nrow(d)) {
    rnk <- rank;
    rnk[is.na(rnk)] <- 0; 
    qn <- abs(qnorm(((1:length(rnk))-0.5)/(2*length(rnk)))); 
    rk <- sign(rnk)*qn[(length(rnk)+1)-abs(rank(abs(rnk)))]
    names(rk) <- names(rnk) <- rownames(d); 
    
    g <- lapply(g[1:2], function(g) g[g %in% colnames(d)]);
    e <- d[, c(g[[1]], g[[2]])];
    e[, g[[1]]] <- apply(e[, g[[1]], drop=FALSE], 2, function(x) x-rowMeans(e[, g[[1]], drop=FALSE])); 
    e[, g[[2]]] <- apply(e[, g[[2]], drop=FALSE], 2, function(x) x-rowMeans(e[, g[[2]], drop=FALSE])+rk);
    e;
  } else {
    d;
  }
}