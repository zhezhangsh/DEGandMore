# EBSeq
# http://bioconductor.org/packages/release/bioc/html/EBSeq.html
# EBTest {EBSeq}
DeEbSeq <- function(mtrx, grps,  paired=FALSE, iter=5, normalization=c('median', 'rank', 'uq', 'none')) {
  require(DEGandMore);
  require(EBSeq);

  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by EBSeq; performing unpaired test instead.\n");
  
  nm <- tolower(normalization)[1];
  
  if (nm == 'none') sz <- rep(1, ncol(mtrx)) else 
    if (nm == 'uq') sz <- QuantileNorm(mtrx, 0.75) else 
      if (nm == 'rank') sz <- RankNorm(mtrx) else 
        sz <- MedianNorm(mtrx); 

  sz <- as.vector(sz); 
  f  <- as.factor(rep(names(grps), sapply(grps, length)));
  eb <- EBTest(Data=mtrx, Conditions = f, sizeFactors = sz, maxround = iter, Print = FALSE); 
  
  pp <- GetPPMat(eb); 
  pv <- pp[, 1][rownames(mtrx)];
  l2 <- pp[, 2][rownames(mtrx)]; 
  
  d  <- eb$DataNorm; 
  m1 <- eb$C1Mean[[1]];
  m2 <- eb$C2Mean[[1]];
  l2 <- log2(pmax(0.5, m2)) - log2(pmax(0.5, m1)); 
  pv <- 1-eb$PPDE;
  qv <- p.adjust(pv, method='BH');
  
  names(l2) <- names(m1); 
  
  id <- rownames(mtrx); 
  m1 <- m1[id];
  m2 <- m2[id];
  l2 <- l2[id];
  pv <- pv[id];
  qv <- qv[id]; 
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s[rownames(mtrx), ], group=grps, eb=eb); 
}