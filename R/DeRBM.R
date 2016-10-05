# RBM
# http://bioconductor.org/packages/3.3/bioc/html/RBM.html
# RBM_T {RBM}

DeRBM <- function(mtrx, grps, paired=FALSE, logged=TRUE, iter=1, ...) {
  
  require(DEGandMore);
  require(RBM);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by RBM; performing unpaired test instead.\n");
  
  f   <- rep(0:1, sapply(grps, length)); 
  res <- RBM_T(mtrx, f, repetition=iter); 
  
  pv <- res$ordfit_pvalue; 
  qv <- p.adjust(pv, method='BH');
  m1 <- res$ordfit_beta0
  m2 <- m1 + res$ordfit_beta1; 
  if (logged) l2 <- m2-m1 else {
    mn <- min(mtrx[mtrx>0]); 
    l2 <- log2(pmax(mn, m2)) - log2(pmax(mn, m1)); 
  };
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  s <- cbind(s, permutation_p = res$permutation_p, bootstrap_p = res$bootstrap_p); 

  list(stat=s, group=grps, RBM=res);
}