# ABSSeq
# https://bioconductor.org/packages/release/bioc/html/ABSSeq.html
# callDEs {ABSSeq}
DeAbsSeq <- function(mtrx, grps, paired=FALSE, norm.method=c("geometric", "quartile", "total", "none"), ...) {
  # norm.method   suggestion: use 'upperquantile' for ChIP-seq data
  
  require(DEGandMore);
  require(ABSSeq);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  nm <- tolower(norm.method)[1];
  if (!(nm %in% c("geometric", "quartile", "total"))) {
    nm <- 'user';
    sz <- rep(1, ncol(mtrx)); 
  } else sz <- 0; 
  f <- factor(rep(names(grps), sapply(grps, length))); 
  
  obj <- ABSDataSet(mtrx, groups = f, normMethod = nm, sizeFactor = sz, paired=paired);
  obj <- ABSSeq(obj); 
  res <- ABSSeq::results(obj)[rownames(mtrx), , drop=FALSE];
  
  m1 <- exp(log(2)*res[, 'Amean']) - 1;
  m2 <- exp(log(2)*res[, 'Bmean']) - 1;
  l2 <- res[, 'Bmean'] - res[, 'Amean']; 
  pv <- res[, 'pvalue']; 
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
  
  if(is.loaded('ABSSeq')) detach("package:ABSSeq", unload=TRUE);
  
  list(stat=s[rownames(mtrx), ], group=grps, ABSSeq=obj);
}
  