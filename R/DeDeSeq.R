DeDeSeq <- function(mtrx, grps, paired=FALSE, ...) {
  require(DEGandMore);
  require(DESeq2);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  n <- sapply(grps, length);
  
  # unpaired test
  if (paired) {
    cond <- DataFrame(cond=factor(rep(names(grps), n)), pair=factor(c(1:n[1], 1:n[2])));
    dds <- DESeqDataSetFromMatrix(mtrx, cond, ~ cond + pair);
  } else {
    cond <- DataFrame(cond=factor(rep(names(grps), n))); 
    dds  <- DESeqDataSetFromMatrix(mtrx, cond, ~ cond);
  }
      
  dds  <- DESeq(dds);
  res  <- results(dds, c('cond', names(grps)[2:1]));
  
  # Summary statistics
  lgfc <- res[, 2];
  pval <- 1-2*abs(pnorm(res[, 4])-0.5); 
  lgfc[is.na(lgfc)] <- 0;
  pval[is.na(pval)] <- 1;  
  
  nm <- sapply(1:ncol(mtrx), function(i) mtrx[, i]/dds@colData@listData$sizeFactor[i]); 
  e1 <- nm[, grps[[1]], drop=FALSE];
  e2 <- nm[, grps[[2]], drop=FALSE];
  m1 <- rowMeans(e1, na.rm=TRUE);
  m2 <- rowMeans(e2, na.rm=TRUE);
  q <- p.adjust(pval, method='BH');
  s <- cbind(m1, m2, m2-m1, lgfc, pval, q);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps)[2:1], collapse='-'), 
                   'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  

  list(stat=s, group=grps, dds=dds);
}