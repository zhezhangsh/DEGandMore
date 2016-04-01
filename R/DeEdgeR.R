DeEdgeR<-function(mtrx, grps, norm.method='TMM', ...) {
  # norm.method   suggestion: use 'upperquantile' for ChIP-seq data
  library(DEGandMore);
  library(edgeR);
  
  e1<-mtrx[, grps[[1]], drop=FALSE];
  e2<-mtrx[, grps[[2]], drop=FALSE];
  ct<-cbind(e1, e2);
  group<-rep(names(grps), sapply(grps, length));
  
  # create DGEList object, normalize data and estimate dispersion 
  dge<-DGEList(counts=ct, group=group);
  dge<-calcNormFactors(dge, method=norm.method);
  dge<-estimateCommonDisp(dge);
  if (ncol(ct)==2) dge@.Data[[3]]<-0.5 else # No replicates
    dge<-estimateTagwiseDisp(dge); 
  
  stat<-as.data.frame(exactTest(dge)[[1]]);
  m1<-rowMeans(e1, na.rm=TRUE);
  m2<-rowMeans(e2, na.rm=TRUE);
  lgfc<-stat[,1];
  lgfc[is.na(lgfc)]<-0;
  p<-stat[,3];
  p[is.na(p)]<-1;
  q<-p.adjust(p, method='BH');
  
  s<-cbind(m1, m2, m2-m1, lgfc, p, q);
  colnames(s)<-c(paste('Mean', names(grps), sep='_'), paste(names(grps), collapse='-'), 'LogFC', 'Pvalue', 'FDR');
  s<-cbind(s, stat[, 2, drop=FALSE]);
  
  list(stat=s[rownames(mtrx), ], group=grps, dge=dge);
}