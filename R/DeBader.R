# BADER
# https://www.bioconductor.org/packages/3.3/bioc/html/BADER.html
# BADER {BADER}
DeBader <- function(mtrx, grps, paired=FALSE, normalized=FALSE, burn = 1000, reps = 10000, ...) {
  require(DEGandMore);
  require(BADER);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by BADER; performing unpaired test instead.\n");
  
  design <- factor(rep(names(grps), sapply(grps, length)));  
  
  res <- BADER(mtrx, design, sizeFactors = !normalized, burn = burn, reps = reps, 
                    printEvery=1000, saveEvery = 1000, mode = 'full'); 
  
  m1 <- res$logMeanA;
  m2 <- res$logMeanB;
  l2 <- log2(exp(res$logFoldChange)); 
  pv <- 1 - res$diffProb; 
  qv <- p.adjust(pv, method='BH');
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps)[2:1], collapse='-'), 
                   'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, bader=res);
}

