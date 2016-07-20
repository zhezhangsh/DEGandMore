# A default run of Rank Product test on 2 groups of samples, and summary of results

###################################################################
# Make output consistent with other methods
DeRankP<-function(mtrx, grps, paired=FALSE, logged=TRUE, nperm=min(100, 2*ncol(mtrx)), ...) {
  require(RankProd);
  require(DEGandMore);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  stat<-DeRankProd(mtrx=mtrx, grps=grps, paired=paired, logged=logged, nperm=nperm, ...);
  
  s<-stat$single.rank[, 4:6];
  if (logged) s<-cbind(s, s[, 3]) else s<-cbind(s, log2(s[,2]/s[,1]));
  s[is.na(s)]<-0;
  s<-cbind(s, stat$single.rank[, c(8, 9, 2, 3, 1)]);
  s[, 'FDR'] <- p.adjust(s[, 'Pvalue'], method='BH'); 
  
  colnames(s)<-c(paste('Mean', names(grps), sep='_'), paste(names(grps), collapse='-'), 'LogFC', 'Pvalue', 'FDR', 'RP1', 'RP2', 'Rank');
  
  list(stat=s[rownames(mtrx), ], group=grps, rp=list(original=stat$rp, summarized=stat[1:2]));
}

DeRankProd<-function(mtrx, grps, paired=FALSE, logged=TRUE, nperm=100, ...) {
  require(RankProd);
  
  nm<-names(grps);
  default.nm<-c('Control', 'Case');
  if (is.null(nm)) nm<-default.nm else nm[is.na(nm) | nm=='']<-default.nm[is.na(nm) | nm==''];
  
  if (paired & length(grps[[1]])==length(grps[[2]])) {
    capture.output(rp<-RP(mtrx[, grps[[1]]]-mtrx[, grps[[2]]], rep(1, length(grps[[1]])), num.perm=nperm, logged=logged, ...))->x;
  } else {
    capture.output(rp<-RP(mtrx[, c(grps[[1]], grps[[2]])], rep(0:1, sapply(grps, length)), num.perm=nperm, logged=logged, ...))->x;
  }

  stat<-SummarizeRP(rp, nm[1], nm[2], save.it=FALSE, single.ranking=TRUE, write.rnk=FALSE, write.excel=FALSE);
  
  means<-sapply(grps, function(c) rowMeans(mtrx[, c, drop=FALSE]));
  colnames(means)<-paste('Mean', nm, sep='_');
  
  cnm<-c(paste("Log2(", nm[2], '/', nm[1], ')', sep=''), "FoldChange", "Pvalue", "FDR");
  stat<-lapply(stat, function(stat) cbind(stat[, !(colnames(stat) %in% cnm)], means, stat[, cnm]));
  
  stat$rp<-rp;
  
  stat;
}

#################################################################################################################################
# Summarize RP results
SummarizeRP<-function(rp, class1, class2, genename=NA, save.it=TRUE, single.ranking=FALSE, write.rnk=FALSE, write.excel=FALSE) {
  
  tbs<-lapply(1:2, function(i) cbind(sapply(rp[4:1], function(rp) rp[[i]]), -1*rp$AveFC, exp(-1*rp$AveFC*log(2))));
  colnames(tbs[[1]])<-colnames(tbs[[2]])<-c('Rank', 'RankProduct', 'Pvalue', 'FDR', paste('Log2(', class2, '/', class1, ')', sep=''), 'FoldChange');
  if (!identical(NA, genename)) rownames(tbs[[1]])<-rownames(tbs[[2]])<-genename;
  names(tbs)[1]<-paste(class1, '<', class2);
  names(tbs)[2]<-paste(class1, '>', class2);
  ps <- rp$pval; 
  fc <- rp$AveFC;
  
  if (single.ranking) { # use a single ranking for both directions of change
    rk<-rank(log2(tbs[[1]][,2])+log2(1/tbs[[2]][,2]));
    p <- ps[, 1];
    p[fc > 0] <- ps[fc > 0, 2]; 
    p[is.na(fc) | is.na(p)] <- 1;
#     p<-rep(1, length(rk)); 
#     ind1<-tbs[[1]][, 5]>0 & !is.na(tbs[[1]][, 5]); 
#     p[ind1]<-tbs[[1]][ind1, 3];
#     ind2<-tbs[[2]][, 5]>0 & !is.na(tbs[[2]][, 5]); 
#     p[ind2]<-tbs[[2]][ind2, 3];
    p<-pmin(1, 2*p);
    a<-cbind(Rank=rk, RP1=tbs[[1]][,2], RP2=tbs[[2]][,2], tbs[[1]][, 5:6], Pvalue=p, FDR=p.adjust(p, method='BH'));
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