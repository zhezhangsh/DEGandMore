DeDeSeq2 <- function(mtrx, grps, paired=FALSE, ...) {
  require(DEGandMore);
  require(DESeq2);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  e1 <- mtrx[, grps[[1]], drop=FALSE];
  e2 <- mtrx[, grps[[2]], drop=FALSE];
  
  colData <- data.frame(row.names=colnames(mtrx), Condition=rep(names(grps), sapply(grps, length)), stringsAsFactors = FALSE);
  
  dds <- DESeqDataSetFromMatrix(countData = cbind(e1, e2), colData = colData, design = ~ Condition);
  ds  <- DESeq2::DESeq(dds, quiet=TRUE, ...);
  res <- data.frame(results(ds))[rownames(mtrx), , drop=FALSE];
  
  m1 <- rowMeans(e1, na.rm=TRUE);
  m2 <- rowMeans(e2, na.rm=TRUE);
  stat <- cbind(m1, m2, m2-m1, res[,2], res[, 5], res[, 6]);
  stat[is.na(stat[,4]), 4] <- 0;
  stat[is.na(stat[,5]), 5] <- 1;
  stat[is.na(stat[,6]), 6] <- 1;
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  stat<-cbind(stat, res[, c(1, 3, 4)]);
  
  list(stat=stat[rownames(mtrx), ], group=grps, deseq=ds);
}