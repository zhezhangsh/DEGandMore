# Internal function that calculates log2(fold change) based on normalized group means and original read count matrix
CalculateCountLog2FC <- function(mean0, mean1, mtrx, grps) {
  
  mn0 <- min(mean0[mean0>0 & !is.na(mean0)]);
  mn1 <- min(mean1[mean1>0 & !is.na(mean1)]);
  
  nm0 <- rownames(mtrx)[mean0==mn0 & !is.na(mean0)];
  nm1 <- rownames(mtrx)[mean1==mn1 & !is.na(mean1)];
  
  adj0 <- 0.5/sum(mtrx[nm0, grps[[1]]]);
  adj1 <- 0.5/sum(mtrx[nm1, grps[[2]]]);
  
  log2(pmax(adj0, mean0) / pmax(adj1, mean1));
}