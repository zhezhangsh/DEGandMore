DeTSPM <- function(mtrx, grps, paired=FALSE) {
  require(DEGandMore);
  require(ShrinkBayes); 

  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  n   <- sapply(grps, length); 
 
  if (paired & n[1]==n[2]) {
    
  } else {
    des <- factor(rep(1:2, n));
    frm <- ~ 1 +  des;
    res <- ShrinkBayesWrap(dat=e[1:100, ], form=frm, ncpus2use = 2);
  }
  ##########################################
  res <- TSPM(mtrx, x1, x0, colSums(mtrx)); 
  ##########################################
  
  l2 <- -res$log.fold.change;
  l2 <- log2(exp(l2)); 
  pv <- res$pvalues;
  qv <- p.adjust(pv, method='BH');
  
  aj <- colSums(mtrx)/mean(colSums(mtrx));
  m1 <- rowMeans(mtrx[, grps[[1]]])/mean(aj[grps[[1]]]);
  m2 <- rowMeans(mtrx[, grps[[2]]])/mean(aj[grps[[2]]]);
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, TSPM=res);
}