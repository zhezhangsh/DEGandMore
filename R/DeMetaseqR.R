DeMetaseqR <- function(mtrx, grps, paired=FALSE) {
  require(DEGandMore);
  require(metaseqR); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by metaseqR; performing unpaired test instead.\n");
  
  g <- lapply(grps, function(i) colnames(mtrx)[i]); 

  # baySeq
  norm1 <- normalize.deseq(mtrx, g);
  p.bayseq <- stat.bayseq(norm1, g)[[1]];
  
  # deseq
  p.deseq <- stat.deseq(norm1, g)[[1]];
  
  # edger
  norm2 <- normalize.edger(mtrx, g);
  p.edger <- stat.edger(norm2, g)[[1]];
  
  # limma
  p.limma <- stat.limma(norm2, g)[[1]];
  
  # nbpseq
  norm3 <- normalize.nbpseq(mtrx, g);
  p.nbpseq <- stat.nbpseq(norm3, g)[[1]]; 
  
  # noiseq
  p.noiseq <- stat.noiseq(norm2, g)[[1]]; 
  
  ps <- cbind(p.bayseq, p.deseq, p.edger, p.limma, p.nbpseq, p.noiseq); 
  pv <- apply(ps, 1, combine.simes); 
  qv <- p.adjust(pv, method='BH');
  nm <- (norm1+norm2+norm3)/3;
  m1 <- rowMeans(nm[, grps[[1]]]); 
  m2 <- rowMeans(nm[, grps[[2]]]);
  l2 <- log2(m2+1) - log2(m1+1); 
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, pvalue=ps);
}