# BGmix
# https://www.bioconductor.org/packages/3.3/bioc/html/BGmix.html
# BGmix {BGmix} 
DeBGmix <- function(mtrx, grps, paired=FALSE, logged=TRUE, ...) {
  require(DEGandMore);
  require(BGmix);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  n        <- as.vector(sapply(grps, length)); 
  
  m1 <- rowMeans(mtrx[, grps[[1]], drop=FALSE], na.rm=TRUE); 
  m2 <- rowMeans(mtrx[, grps[[2]], drop=FALSE], na.rm=TRUE); 
  v1 <- apply(mtrx[, grps[[1]], drop=FALSE], 1, function(x) var(x, na.rm=TRUE)); 
  v2 <- apply(mtrx[, grps[[2]], drop=FALSE], 1, function(x) var(x, na.rm=TRUE)); 
  
  if (paired & n[1]==n[2]) {
    m0 <- matrix(m2-m1, nc=1);
    d0 <- mtrx[, grps[[1]], drop=FALSE]-mtrx[, grps[[2]], drop=FALSE]; 
    v0 <- matrix(apply(d0, 1, function(x) var(x, na.rm=TRUE)), nc=1); 
    xx <- matrix(1, nr=1); 
    dir <- BGmix(m0, v0, n[1], xx=xx, neffects=1, jstar = 0, ntau = 1, indtau = 0, niter = 100, nburn = 100); 
    
    par <- ccParams(dir);  
    res <- ccTrace(dir);
    
    var1 <- res$sig2;
    vart <- var1/n[1]; 
    tp <- rowMeans((m2-m1)/sqrt(vart));
    pv <- pnorm(tp); 
    pv[pv>.5] <- 1-pv[pv>.5];
    pv <- 2*pv;
  } else {
    mm <- cbind(m1, m2);
    ss <- cbind(v1, v2);
    dimnames(mm) <- dimnames(ss) <- list(NULL, NULL);
    
    dir <- BGmix(mm, ss, n, niter = 100, nburn = 100); 
    
    par <- ccParams(dir); 
    res <- ccTrace(dir); 
    
    var1 <- res$sig2[1, ,];
    var2 <- res$sig2[2, ,];
    vart <- var1/n[1] + var2/n[2]; 
    tp <- rowMeans((m2-m1)/sqrt(vart));
    pv <- 2*pnorm(-abs(tp));
  }
  
  qv <- p.adjust(pv, method='BH');
  if (logged) l2 <- m2-m1 else {
    mn <- min(mtrx[mtrx>0]); 
    l2 <- log2(pmax(mn, m2)) - log2(pmax(mn, m1)); 
  };
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  unlink(dir, recursive = TRUE); 
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, bgmix=list(parameter=par, trace=res));
}


