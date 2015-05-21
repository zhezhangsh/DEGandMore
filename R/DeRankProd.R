# A default run of Rank Product test on 2 groups of samples, and summary of results
DeRankProd<-function(mtrx, clss, nperm=100) {
  library(RankProd);
  
  nm<-names(clss);
  default.nm<-c('Control', 'Case');
  if (is.null(nm)) nm<-default.nm else nm[is.na(nm) | nm=='']<-default.nm[is.na(nm) | nm==''];
  
  cl<-rep(0:1, sapply(clss, length));
  
  rp<-RP(mtrx[, c(clss[[1]], clss[[2]]), drop=FALSE], cl, num.perm=nperm);
  
  stat<-SummarizeRP(rp, nm[1], nm[2], save.it=FALSE, single.ranking=TRUE, write.rnk=FALSE, write.excel=FALSE);
  
  means<-sapply(clss, function(c) rowMeans(mtrx[, c, drop=FALSE]));
  colnames(means)<-paste('Mean', nm, sep='_');
  
  cnm<-c(paste("Log2(", nm[2], '/', nm[1], ')', sep=''), "FoldChange", "PValue", "FDR")
  stat<-lapply(stat, function(stat) cbind(stat[, !(colnames(stat) %in% cnm)], means, stat[, cnm]));
  
  stat;
}

#################################################################################################################################
# Summarize RP results
SummarizeRP<-function(rp, class1, class2, genename=NA, save.it=TRUE, single.ranking=FALSE, write.rnk=FALSE, write.excel=FALSE) {
  
  tbs<-lapply(1:2, function(i) cbind(sapply(rp[4:1], function(rp) rp[[i]]), -1*rp$AveFC, exp(-1*rp$AveFC*log(2))));
  colnames(tbs[[1]])<-colnames(tbs[[2]])<-c('Rank', 'RankProduct', 'PValue', 'FDR', paste('Log2(', class2, '/', class1, ')', sep=''), 'FoldChange');
  if (!identical(NA, genename)) rownames(tbs[[1]])<-rownames(tbs[[2]])<-genename;
  names(tbs)[1]<-paste(class1, '<', class2);
  names(tbs)[2]<-paste(class1, '>', class2);
  
  if (single.ranking) { # use a single ranking for both directions of change
    rk<-rank(log2(tbs[[1]][,2])+log2(1/tbs[[2]][,2]));
    p<-(tbs[[1]][,3]+(1-tbs[[2]][,3]))/2;
    p<-2*pmin(p, 1-p);
    p[p<0]<-0;
    a<-cbind(Rank=rk, RP1=tbs[[1]][,2], RP2=tbs[[2]][,2], tbs[[1]][, 5:6], PValue=p, FDR=p.adjust(p, method='BH'));
    if (!identical(NA, genename)) rownames(a)<-genename;
    tbs$single.rank<-a;
  }
  
  #if (write.excel) {
  #	xls<-tbs[1:2]
  #	if (single.ranking) xls[['Single rank']]<-cbind(Gene=names(tbs[[3]]), Rank=as.vector(tbs[[3]]));
  #	Excel(xls, paste(class1, 'vs', class2));
  #}
  
  if (write.rnk) {
    if (single.ranking) {
      rnk<-cbind(genename, log2(1/tbs[[1]][,2])+log2(tbs[[2]][,2]));
      rnk<-rnk[!is.na(genename), ];
      write.table(rnk, paste(class1, '-vs-', class2, '.rnk', sep=''), sep='\t', row=FALSE, col=FALSE, qu=FALSE);
    } else {
      write.table(tbs[[1]], paste(class1, '<', class2, '.rnk', sep=''), sep='\t', row=FALSE, col=FALSE, qu=FALSE);
      write.table(tbs[[2]], paste(class1, '>', class2, '.rnk', sep=''), sep='\t', row=FALSE, col=FALSE, qu=FALSE);
    }
  }		
  
  # FDR cannot be greater than 1
  for (i in 1:length(tbs)) tbs[[i]][, 'FDR']<-pmin(1, tbs[[i]][, 'FDR']);
  
  if (save.it) save(tbs, file=paste(class1, '-vs-', class2, '.rdata', sep=''));
  tbs;
}