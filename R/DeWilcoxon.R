# Paired/Unpaired Wilcoxon/Mann-Whitney test
DeWilcoxon <- function(mtrx, grps, paired=FALSE, logged=TRUE) {
  # mtrx    A numeric matrix of gene expression data. Rows are unique genes and columns include 2 groups of samples to be compared
  # grps    A 2-vector list, each vector has the column indexes (numeric vectors) or column names (character vectors) of a group; vectors are named by group names  
  # paired  Whether it's a paired test; if TRUE the two groups must have the same number of samples and the samples must be ordered to match each other
  # logged  Whether data is log2-transformed. It will affect how Log(FoldChange), but not the p values, will be calculated. It won't affect methods using RNA-seq read count data either.
  
  require(DEGandMore);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  # data subsets of the sample groups
  d1 <- mtrx[, grps[[1]], drop=FALSE];
  d2 <- mtrx[, grps[[2]], drop=FALSE];
  
  # sample sizes
  n1 <- ncol(d1);
  n2 <- ncol(d2);
  
  # group means
  m1 <- rowMeans(d1, na.rm=TRUE);
  m2 <- rowMeans(d2, na.rm=TRUE);
  diff <- m2-m1;
  
  if (logged) lgfc <- diff else lgfc <- log2(m2/pmax(m1, min(m1[m1>0])));
  
  p <- sapply(1:nrow(mtrx), function(i) {
    if (!paired | n1!=n2) wilcox.test(d1[i, ], d2[i, ])$p.value else wilcox.test(d1[i, ], d2[i, ], paired=TRUE)$p.value
  }); 
  p[is.na(p)] <- 1;
  stat <- cbind(m1, m2, diff, lgfc, p, p.adjust(p, method='BH'));
  colnames(stat) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(stat) <- rownames(mtrx);

  list(stat=stat[rownames(mtrx), ], group=grps);
}