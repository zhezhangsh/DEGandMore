DeVoomLimma<-function(mtrx, grps, paired=FALSE, plot=FALSE, ...) {
  require(DEGandMore);
  require("edgeR");
  require("limma");

  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  res  <- VoomLimma(mtrx=mtrx, grps=grps, paired=paired, plot=plot, ...);
  
  e1   <- mtrx[, grps[[1]], drop=FALSE];
  e2   <- mtrx[, grps[[2]], drop=FALSE];
  m1   <- rowMeans(e1, na.rm=TRUE);
  m2   <- rowMeans(e2, na.rm=TRUE);

  lgfc <- res[[1]][, 1];
  p    <- res[[1]][, 2];
  q    <- res[[2]][, 'adj.P.Val'];
  lgfc[is.na(lgfc)] <- 0;
  p[is.na(p)]       <- 1;
  q[is.na(q)]       <- 1;
  
  s <- cbind(cbind(m1, m2, m2-m1)[rownames(res[[1]]), ], lgfc, p, q);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps), collapse='-'), 'LogFC', 'Pvalue', 'FDR');
  s <- cbind(s, res[[2]][, c(2, 3, 6)]);
  
  list(stat=s[rownames(mtrx), ], group=grps, voom=res);
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
    x      <- factor(paste('Pair', rep(1:n, 2), sep='_'));
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
  fit2 <- eBayes(fit2)
  
  topTable(fit2, coef='yGroup_2', number=nrow(ct)); 
}
