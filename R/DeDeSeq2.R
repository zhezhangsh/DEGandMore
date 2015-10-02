DeDeSeq2<-function(mtrx, grps, ...) {
  library(DEGandMore);
  library(DESeq2);
  
  e1<-mtrx[, grps[[1]], drop=FALSE];
  e2<-mtrx[, grps[[2]], drop=FALSE];
  
  colData<-data.frame(row.names=unlist(grps, use.names=FALSE), stringsAsFactors = FALSE, Condition=rep(names(grps), sapply(grps, length)));
  
  dds <- DESeqDataSetFromMatrix(countData = cbind(e1, e2), colData = colData, design = ~ Condition);
  ds <- DESeq(dds, quiet=TRUE, ...);
  res<-data.frame(results(ds));
  
  m1<-rowMeans(e1, na.rm=TRUE);
  m2<-rowMeans(e2, na.rm=TRUE);
  stat<-cbind(e1, e2, e2-e1, res[,2], res[, 5], res[, 6]);
  stat[is.na(stat[,5]), 5]<-1;
  stat[is.na(stat[,6]), 6]<-1;
  colnames(stat)<-c(paste('Mean', names(grps), sep='_'), paste(names(grps), collapse='-'), 'LogFC', 'Pvalue', 'FDR');
  stat<-cbind(stat, res[, c(1, 3, 4)]);
  
  list(stat=stat[rownames(mtrx), ], group=grps, deseq=ds);
}