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
  libadjust = apply(mtrx, 2, function(x) {
    lognz = log2(x[x != 0] + 1)
    q3 = quantile(lognz, 0.75)
    sum(lognz[lognz < q3])
  }); 
  libadjust <- libadjust/mean(libadjust); 
  expr <- log2(mtrx+1);
  expr <- sapply(1:ncol(expr), function(i) pmax(log2(4/3), expr[, i])/libadjust[i]); 
  dimnames(expr) <- dimnames(mtrx); 
  
  m1 <- rowMeans(expr[, grps[[1]]]); 
  m2 <- rowMeans(expr[, grps[[2]]]);
  l2 <- m2-m1;
  m1 <- 2^m1;
  m2 <- 2^m2;
  
  if (paired & n[1]==n[2]) {
    pair <- factor(rep(1:n[1], 2)); 
    comp <- factor(rep(names(grps), sapply(grps, length))); 
    phen <- data.frame(group=comp, pair=pair);
    mod0 <- model.matrix(~ phen$pair);
    mod1 <- model.matrix(~ phen$pair + phen$group);
    stat <- stattest(gowntable=expr, pData=phen, mod=mod1, mod0=mod0, feature='gene', log=TRUE);
  } else {
    stat <- stattest(gowntable=expr, pData=phen, covariate='group', feature='gene', log=TRUE);
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
  
  list(stat=s, group=grps);
}