DeNBPSeq <- function(mtrx, grps, paired=FALSE) {
  require(DEGandMore);
  require(NBPSeq); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by NBPSeq; performing unpaired test instead.\n");
  
  n <- sapply(grps, length); 
  
  nm = estimate.norm.factors(mtrx);
  
  ###########################################################################
  res <- nbp.test(mtrx, rep(1:2, n), 1, 2, nm); 
  ###########################################################################
  
  sm <- sum(rowMeans(mtrx)); 
  
  pv <- res$p.values; 
  qv <- p.adjust(pv, method='BH');
  l2 <- res$log.fc; 
  m1 <- sm * res$expression.levels[, 1]
  m2 <- sm * res$expression.levels[, 2]
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, NBPSeq=res);
}