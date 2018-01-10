DeROTS <- function(mtrx, grps, paired=FALSE, logged=TRUE, B=1000) {
  
  require(DEGandMore);
  require(ROTS); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];

  g <- c(rep(0, length(grps[[1]])), rep(1, length(grps[[2]])));  
  k <- min(nrow(mtrx)/2, min(5000, 1000*ceiling(ncol(mtrx)/50000*nrow(mtrx)))); 
  res <- ROTS(mtrx, g, B=B, K=k, log=logged);
  
  m1 <- rowMeans(mtrx[, grps[[1]], drop=FALSE]);
  m2 <- rowMeans(mtrx[, grps[[2]], drop=FALSE]);
  fc <- res$logfc;
  pv <- res$pvalue;
  qv <- res$FDR;
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  fc[is.na(fc)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, fc, pv, qv, res$d);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR', 'Stat_d');
  rownames(s) <- rownames(mtrx); 

  list(stat=s, group=grps, rots=res[c('B', 'a1', 'a2', 'k', 'R', 'Z')]);
}