DeBallgown <- function(mtrx, grps, paired=FALSE) {
  
  require(DEGandMore);
  require(ballgown); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  n <- sapply(grps, length);
  
  # Unpair model
  comp <- factor(rep(names(grps), sapply(grps, length))); 
  phen <- data.frame(group=comp);

  # Calculate log2-FC using the upper-quanitle method. Detail at stattest {ballgown} help page, parameter libadjust
  lg <- log2(mtrx+1); 
  l2 <- log2(stattest(gowntable=lg, pData=phen, covariate='group', feature='gene', getFC=TRUE, log=TRUE)[, 3]);
  m1 <- rowMeans(mtrx[, grps[[1]]]); 
  m2 <- m1*2^l2;
  
  if (paired & n[1]==n[2]) {
    pair <- factor(rep(1:n[1], 2)); 
    comp <- factor(rep(names(grps), sapply(grps, length))); 
    phen <- data.frame(group=comp, pair=pair);
    mod0 <- model.matrix(~ phen$pair);
    mod1 <- model.matrix(~ phen$pair + phen$group);
    stat <- stattest(gowntable=mtrx, pData=phen, mod=mod1, mod0=mod0, feature='gene', log=FALSE);
  } else {
    stat <- stattest(gowntable=mtrx, pData=phen, covariate='group', feature='gene', log=FALSE);
  }

  pv <- stat[, 3]; 
  qv <- p.adjust(pv, method='BH');

  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, limma=fit);
}