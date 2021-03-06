# PoissonSeq
# https://cran.r-project.org/web/packages/PoissonSeq/index.html
# PS.Main {PoissonSeq}
DePoissonSeq <- function(mtrx, grps, paired = FALSE) {
  
  require(DEGandMore);
  require(PoissonSeq);

  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  dat <- list(n = mtrx, y = rep(1:2, sapply(grps, length)), type = 'twoclass', pair = paired);
  par <- list(trans = TRUE, npermu = min(1000, 20*ncol(mtrx)), ct.sum=0, ct.mean=0); 
  res <- PS.Main(dat, para=par); 
   
  norm <- NormTMM(mtrx); # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2864565/ https://www.ncbi.nlm.nih.gov/pubmed/20196867/
  
  m1 <- rowMeans(norm[, grps[[1]], drop=FALSE]);
  m2 <- rowMeans(norm[, grps[[2]], drop=FALSE]);
  l2 <- res[rownames(mtrx), 'log.fc']; 
  pv <- res[rownames(mtrx), 'pval']; 
  qv <- p.adjust(pv, method='BH');
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s[rownames(mtrx), ], group=grps, poissonseq=list(dat=dat, para=par)); 
}