# Significance analysis of sequencing data
# SAMseq {samr}
DeSamSeq <- function(mtrx, grps, paired=FALSE, normalization=c('TMM', 'RLE', 'DESeq', 'Median', 'UQ', 'TC', 'QQ'), 
                     nperms=max(100, 10*ncol(mtrx)), nresamp=max(20, 4*ncol(mtrx)), ...) {
  # norm.method   suggestion: use 'upperquantile' for ChIP-seq data
  require(DEGandMore);
  require(samr);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  
  norm <- tolower(normalization)[1]; 
  if (normalization=='tmm')    mtrx <- NormTMM(mtrx);
  if (normalization=='rle')    mtrx <- NormRLE(mtrx);
  if (normalization=='deseq')  mtrx <- NormDESeq(mtrx);
  if (normalization=='median') mtrx <- NormMedian(mtrx);
  if (normalization=='uq')     mtrx <- NormUpperQuantile(mtrx);
  if (normalization=='tc')     mtrx <- NormTotalCount(mtrx);
  if (normalization=='qq')     mtrx <- NormQQ(mtrx);
  mtrx <- round(mtrx); 
  mtrx <- mtrx[, c(grps[[1]], grps[[2]]), drop=FALSE]; 
  
  n    <- sapply(grps, length); 
  
  capture.output(
    if (paired & n[1]==n[2]) {
      y <- c(-(1:n[1]), 1:n[1]);
      sam <- SAMseq(mtrx, y, resp.type = 'Two class paired', nperms = nperms, nresamp = nresamp)$samr.obj;
    } else {
      y <- c(rep(1, n[1]), rep(2, n[2])); 
      sam <- SAMseq(mtrx, y, resp.type = 'Two class unpaired', nperms = nperms, nresamp = nresamp)$samr.obj;
  }) -> x;
  
  d  <- list(x=mtrx, y=y, geneid=rownames(mtrx), genename=rownames(mtrx), logged2=FALSE); 
  p  <- samr.pvalues.from.perms(sam$tt, sam$ttstar);
  m1 <- rowMeans(mtrx[, grps[[1]], drop=FALSE], na.rm=TRUE);
  m2 <- rowMeans(mtrx[, grps[[2]], drop=FALSE], na.rm=TRUE);
  q  <- p.adjust(p, method='BH');
  fc <- log2(pmax(0.5, m2)) - log2(pmax(0.5, m1)); 
  s  <- cbind(m1, m2, m2-m1, fc, p, q);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps)[2:1], collapse='-'), 
                   'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s[rownames(mtrx), ], group=grps, cd=cd);
}