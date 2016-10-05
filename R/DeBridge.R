# bridge
# http://bioconductor.org/packages/3.3/bioc/html/bridge.html
# DeBridge
DeBridge <- function(mtrx, grps, paired=FALSE, logged=TRUE, ...) {
  require(DEGandMore);
  require(bridge);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by this method; performing unpaired test instead.\n");
  
  d1  <- mtrx[, grps[[1]], drop=FALSE];
  d2  <- mtrx[, grps[[2]], drop=FALSE];
  
  res <- bridge.2samples(d1, d2, log=logged, all.out=FALSE);
  
  pv  <- 1-res$post.p; 
  qv <- p.adjust(pv, method='BH');
  m1 <- rowMeans(d1, na.rm=TRUE);
  m2 <- rowMeans(d2, na.rm=TRUE);
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
  
  list(stat=s, group=grps, bridge=res);
}