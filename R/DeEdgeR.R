DeEdgeR <- function(mtrx, grps, paired=FALSE, norm.method='TMM', ...) {
  # norm.method   suggestion: use 'upperquantile' for ChIP-seq data
  require(DEGandMore);
  require(edgeR);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  n <- sapply(grps, length);
  
  # create DGEList object, normalize data and estimate dispersion 
  group <- rep(names(grps), sapply(grps, length));
  dge   <- DGEList(counts=mtrx, group=group);
  dge   <- calcNormFactors(dge, method=norm.method);
  
  if (paired & n[1]==n[2]) {
    n <- length(grps[[1]]); 
    x <- factor(paste('Pair', rep(1:n, 2), sep='_'));
    y <- factor(paste('Group', rep(1:2, each=n), sep='_'));
    
    design <- model.matrix(~x+y);
    dge    <- estimateGLMCommonDisp(dge, design, verbose=FALSE);
    dge    <- estimateGLMTrendedDisp(dge, design);
    dge    <- estimateGLMTagwiseDisp(dge, design);
    
    fit  <- glmFit(dge, design);
    lrt  <- glmLRT(fit);
    stat <- as.data.frame(topTags(lrt, n=nrow(mtrx))); 
  } else {
    dge  <- estimateCommonDisp(dge);
    if (ncol(mtrx)==2) dge@.Data[[3]] <- 0.5 else # No replicates
      dge <- estimateTagwiseDisp(dge); 
    stat <- as.data.frame(exactTest(dge, pair=names(grps))[[1]]);
  }
  stat <- stat[rownames(mtrx), ]; 
  
  nm <- NormTMM(mtrx); 
  m1 <- rowMeans(nm[, grps[[1]], drop=FALSE], na.rm=TRUE);
  m2 <- rowMeans(nm[, grps[[2]], drop=FALSE], na.rm=TRUE);
  l2 <- stat[, 'logFC'];
  p  <- stat[, 'PValue'];
  q  <- p.adjust(p, method='BH');
  
  l2[is.na(l2)] <- 0;
  p[is.na(p)]   <- 1;
  q[is.na(q)]   <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, p, q, stat[, 2]);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR', colnames(stat)[2]);
  rownames(s) <- rownames(mtrx); 

  list(stat=s[rownames(mtrx), ], group=grps, dge=dge);
}
