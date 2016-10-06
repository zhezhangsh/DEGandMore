# http://bioinformatics.oxfordjournals.org/content/early/2015/04/21/bioinformatics.btv209
# UCexactTest {edgeRun}
DeEdgeRun <- function(mtrx, grps, paired=FALSE, norm.method='TMM', iterations=50000, ...) {
  # norm.method   suggestion: use 'upperquantile' for ChIP-seq data
  require(DEGandMore);
  require(edgeR);
  require(edgeRun);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  # create DGEList object, normalize data and estimate dispersion 
  group <- rep(names(grps), sapply(grps, length), lib.size=colSums(mtrx)); 
  dge   <- DGEList(counts=mtrx, group=group);
  dge   <- calcNormFactors(dge, method=norm.method);

  if (paired & length(grps[[1]])==length(grps[[2]])) {
    n <- length(grps[[1]]); 
    x <- factor(paste('Pair', rep(1:n, 2), sep='_'));
    y <- factor(paste('Group', rep(1:2, each=n), sep='_'));
    
    design <- model.matrix(~x+y);
    dge    <- estimateGLMCommonDisp(dge, design, verbose=FALSE);
    dge    <- estimateGLMTrendedDisp(dge, design);
    dge    <- estimateGLMTagwiseDisp(dge, design);
  } else {
    dge  <- estimateCommonDisp(dge);
    if (ncol(mtrx)==2) dge@.Data[[3]] <- 0.5 else # No replicates
      dge <- estimateTagwiseDisp(dge); 
  }
  
  fit <- UCexactTest(dge, upper=iterations);  
  res <- topTags(fit, n=nrow(mtrx))[[1]][rownames(mtrx), , drop=FALSE]; 
  
  nm <- dge@.Data[[4]]; 
  m1 <- rowMeans(nm[, grps[[1]], drop=FALSE], na.rm=TRUE);
  m2 <- rowMeans(nm[, grps[[2]], drop=FALSE], na.rm=TRUE);
  lg <- res[, 'logFC'];
  p  <- res[, 'PValue'];
  q  <- p.adjust(p, method='BH');
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  lg[is.na(lg)] <- 0;
  p[is.na(p)]   <- 1;
  q[is.na(q)]   <- 1;
  
  s <- cbind(m1, m2, m2-m1, lg, p, q, res[, 2]);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR', colnames(res)[2]);
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s[rownames(mtrx), ], group=grps, dge=dge);
}
