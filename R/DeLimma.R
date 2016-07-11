DeLimma <- function(mtrx, grps, paired=FALSE, logged=TRUE) {
  
  require(DEGandMore);
  require(limma); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired & n[1]==n[2]) {
    pair   <- factor(rep(1:2, each=n[1])); 
    comp   <- factor(rep(names(grps), sapply(grps, length))); 
    design <- model.matrix(~pair+comp);
    fit    <- lmFit(mtrx, design);
    fit    <- eBayes(fit); 
    res    <- topTable(fit, coef=2, number=nrow(mtrx));
    res    <- res[rownames(mtrx), , drop=FALSE]; 
  } else {
    design <- cbind(Ctrl=rep(1, ncol(mtrx)), Comp=rep(1:2, sapply(grps, length))); 
    fit    <- lmFit(mtrx, design); 
    fit    <- eBayes(fit); 
    res    <- topTable(fit, coef=2, number=nrow(mtrx)); 
    res    <- res[rownames(mtrx), , drop=FALSE]; 
  }
  
  m1 <- rowMeans(mtrx[, grps[[1]], drop=FALSE]);
  m2 <- rowMeans(mtrx[, grps[[2]], drop=FALSE]);
  pv <- res[, 4]; 
  qv <- p.adjust(pv, method='BH');
  if (logged) l2 <- m2-m1 else l2 <- log2(m2/pmax(m1, min(m1[m1>0])/2));
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps)[2:1], collapse='-'), 
                   'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, limma=fit);
}