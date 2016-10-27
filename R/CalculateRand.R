CalculateRand <- function(mtrx, grps, correct=TRUE) {
  require(DEGandMore);
  require(flexclust); 
  
  grps <- lapply(grps, function(x) x[x %in% colnames(mtrx)]); 
  mtrx <- mtrx[, c(grps[[1]], grps[[2]])]; 
  
  g0 <- rep(0:1, sapply(grps, length));
  hc <- hclust(as.dist(1-cor(mtrx, use='pair'))); 
  k2 <- cutree(hc, k=2); 
  
  randIndex(xtabs(~g0+k2), correct = correct); 
}