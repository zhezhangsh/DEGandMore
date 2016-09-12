# tweeDEseq
# http://bioconductor.org/packages/release/bioc/html/tweeDEseq.html
# tweeDE {tweeDEseq}
DeTweeDeSeq <- function(mtrx, grps, paired=FALSE, cl=4, ...) {
  
  require(DEGandMore);
  require(tweeDEseq);
  require(multicore); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by tweeDE; performing unpaired test instead.\n");
  
  f   <- rep(names(grps), sapply(grps, length)); 
  ct  <- normalizeCounts(mtrx); 
  res <- tweeDE(ct, group = f, mc.cores = max(1, cl)); 
  
  pv <- res[, 6];
  m1 <- res[, 2];
  m2 <- res[, 3];
  l2 <- res[, 4];
  qv <- p.adjust(pv, method='BH');
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps)[2:1], collapse='-'), 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s[rownames(mtrx), ], group=grps, tweeDEseq=res);
}