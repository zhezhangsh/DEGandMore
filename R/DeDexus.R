DeDexus <- function(mtrx, grps, paired=FALSE, norm='RLE') {
  require(DEGandMore);
  require(dexus); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by dexus; performing unpaired test instead.\n");
  
  n <- sapply(grps, length);
  l <- as.factor(rep(names(grps), sapply(grps, length))); 
  
  ###########################################################################
  res <- dexus(mtrx, label = l, normalization = norm, resultObject = 'list'); 
  ###########################################################################
  
  pv <- res$pval; 
  qv <- p.adjust(pv, method='BH');
  l2 <- res$logfc[, 2]/log(2);
  m1 <- res$means[, 1]; 
  m2 <- res$means[, 2]; 
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, dexus=res);
}