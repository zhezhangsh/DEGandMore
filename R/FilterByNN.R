# Filter rows of a data matrix by their distance to nearest neighbor
# Distance to all nearest neighbors must be significantly smaller than the background by given number of SD
# Background is established by randomly selecting pairs of row and calculate their distance
FilterByNN <- function(d, k=9, background.size=10000, background.nsd=abs(qnorm(pmin(0.25, 4*k/nrow(d))))) {
  # d     Data matrix with multiple rows
  # k     Number of nearest neighbors
  # background.size   Number of random pairs of rows to draw from data matrix; 
  # background.nsd    Number of standard deviation from mean distance of random pairs
  
  require(FNN);
  
  k <- min(k, nrow(d));
  
  nr <- nrow(d); 
  bg <- sapply(1:background.size, function(i) {
    pr <- d[sample(1:nr, 2), ];
    sqrt(sum((pr[1, ]-pr[2, ])^2, na.rm=TRUE));
  }); 
  cf <- mean(bg) - background.nsd*sd(bg); 
    
  nn <- get.knn(d, k=k); 
  dt <- nn[[2]]; 
  dt <- dt[, ncol(dt)]; 
  
  tb <- cbind(dt, dt < cf); 
  colnames(tb) <- c('Distance', 'Threshold');
  rownames(tb) <- rownames(d); 
  
  tb;
}