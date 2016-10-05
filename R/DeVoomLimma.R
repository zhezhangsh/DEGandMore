DeVoomLimma<-function(mtrx, grps, paired=FALSE, plot=FALSE, ...) {
  require(DEGandMore);
  require("edgeR");
  require("limma");

  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  res  <- VoomLimma(mtrx=mtrx, grps=grps, paired=paired, plot=plot, ...);

  stat <- res[[1]][rownames(mtrx), , drop=FALSE]; 
  norm <- res[[2]][rownames(mtrx), , drop=FALSE];
  
  adj <- mean(mtrx, na.rm=TRUE)/mean(2^norm, na.rm=TRUE);
  
  m1 <- 2^rowMeans(norm[, grps[[1]], drop=FALSE], na.rm=TRUE);
  m2 <- 2^rowMeans(norm[, grps[[2]], drop=FALSE], na.rm=TRUE);
  m1 <- pmax(0, adj*m1-0.5);
  m2 <- pmax(0, adj*m2-0.5);
  
  l2 <- stat[, 1];
  p  <- stat[, 4];
  q  <- p.adjust(p, method='BH');
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  p[is.na(p)]   <- 1;
  q[is.na(q)]   <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, p, q);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, voom=res);
}


VoomLimma <- function(mtrx, grps, paired=FALSE, plot=FALSE, ...) {
  #Call libraries
  require("edgeR");
  require("limma");
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  e1 <- mtrx[, grps[[1]], drop=FALSE];
  e2 <- mtrx[, grps[[2]], drop=FALSE];
  ct <- cbind(e1, e2);
  group <- rep(names(grps), sapply(grps, length));
  
  if (paired) {
    x      <- factor(paste('Pair', rep(1:length(grps[[1]]), 2), sep='_'));
    y      <- factor(paste('Group', rep(1:2, each=length(grps[[1]])), sep='_'));
    design <- model.matrix(~x+y);
  } else {
    y      <- factor(paste('Group', rep(c(1, 2), sapply(grps, length)), sep='_'));
    design <- model.matrix(~y);
  }
  
  dge  <- DGEList(counts=ct, group=group);
  dge  <- calcNormFactors(dge);
  v    <- voom(dge, design, plot=plot)$E;
  fit  <- lmFit(v, design);
  fit2 <- eBayes(fit)
  
  list(stat=topTable(fit2, coef='yGroup_2', number=nrow(ct)), normalized=v); 
}
