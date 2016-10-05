# sSeq
# http://bioconductor.org/packages/3.3/bioc/html/sSeq.html
# nbTestSH {sSeq}
DeSSeq <- function(mtrx, grps, paired=FALSE, ...) {
  
  require(DEGandMore);
  require(sSeq);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  mtrx <- round(mtrx); 
  n    <- sapply(grps, length); 
  cond <- rep(names(grps), n); 
  
  if (paired) {
    warning("Paired test of sSeq is not properly implemented: direction of changes ignored. Running unpaired test instead.\n");
    paired <- FALSE;
  }
  
  if (paired & n[1]==n[2]) {
    des <- data.frame(subjects = rep(1:n[1], 2));
    res <- nbTestSH(mtrx, cond, names(grps)[1], names(grps)[2], coLevels = des, keepLevelsConsistant = TRUE,
                    pairedDesign=TRUE, pairedDesign.dispMethod="pooled"); 
  } else {
    res <- nbTestSH(mtrx, cond, names(grps)[1], names(grps)[2]); 
  }
  
  pv <- res[, 'pval']; 
  m1 <- res[, 2];
  m2 <- res[, 3];
  l2 <- -res[, 4]; 
  qv <- p.adjust(pv, method='BH');
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s[rownames(mtrx), ], group=grps, sseq=res);
}
