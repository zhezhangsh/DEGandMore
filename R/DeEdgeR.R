DeEdgeR <- function(mtrx, grps, paired=FALSE, norm.method='TMM', ...) {
  # norm.method   suggestion: use 'upperquantile' for ChIP-seq data
  require(DEGandMore);
  require(edgeR);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  e1 <- mtrx[, grps[[1]], drop=FALSE];
  e2 <- mtrx[, grps[[2]], drop=FALSE];
  ct <- cbind(e1, e2);
  group <- rep(names(grps), sapply(grps, length));
  
  # create DGEList object, normalize data and estimate dispersion 
  dge <- DGEList(counts=ct, group=group);
  dge <- calcNormFactors(dge, method=norm.method);
  
  if (paired) {
    n <- length(grps[[1]]); 
    x <- factor(paste('Pair', rep(1:n, 2), sep='_'));
    y <- factor(paste('Group', rep(1:2, each=n), sep='_'));
    
    design <- model.matrix(~x+y);
    dge <- estimateGLMCommonDisp(dge, design, verbose=FALSE);
    dge <- estimateGLMTrendedDisp(dge, design);
    dge <- estimateGLMTagwiseDisp(dge, design);
    
    fit <- glmFit(dge, design);
    lrt <- glmLRT(fit);
    stat <- as.data.frame(topTags(lrt, n=nrow(mtrx))); 
  } else {
    dge <- estimateCommonDisp(dge);
    if (ncol(ct)==2) dge@.Data[[3]] <- 0.5 else # No replicates
      dge <- estimateTagwiseDisp(dge); 
    stat <- as.data.frame(exactTest(dge)[[1]]);
  }
  
  stat <- stat[rownames(ct), ]; 
  m1 <- rowMeans(e1, na.rm=TRUE);
  m2 <- rowMeans(e2, na.rm=TRUE);
  lgfc <- stat[, 'logFC'];
  lgfc[is.na(lgfc)] <- 0;
  p <- stat[, 'PValue'];
  p[is.na(p)] <- 1;
  q <- p.adjust(p, method='BH');
  
  s <- cbind(m1, m2, m2-m1, lgfc, p, q);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps)[2:1], collapse='-'), 'LogFC', 'Pvalue', 'FDR');
  s <- cbind(s, stat[, 2, drop=FALSE]);
  
  list(stat=s[rownames(mtrx), ], group=grps, dge=dge);
}
