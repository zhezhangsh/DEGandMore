DeGPseq <- function(mtrx, grps, paired=FALSE) {
  require(DEGandMore);
  require(GPseq); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by GPSeq; performing unpaired test instead.\n");
  
  u1 <- summary(rowMeans(mtrx[rowSums(mtrx)>0, grps[[1]]]))[[5]];
  u2 <- summary(rowMeans(mtrx[rowSums(mtrx)>0, grps[[2]]]))[[5]];
  w  <- as.vector(u2)/as.vector(u1); 
  
  chi <- apply(mtrx, 1, function(d) {
    if (max(d) > 40000) d <- d / (max(d)/40000); 
    dx <- d[grps[[1]]];
    dy <- d[grps[[2]]];
    ox <- generalized_poisson_likelihood(dx);
    oy <- generalized_poisson_likelihood(dy);
    likelihood_ratio_tissue_generalized_poisson(dx, ox$lambda, ox$theta, dy, oy$lambda, oy$theta, w=w)[[2]]; 
  }); 
  
  pv <- pchisq(chi, 1, lower.tail = FALSE, log.p = TRUE); 
  pv <- exp(pv); 
  qv <- p.adjust(pv, method='BH');

  m1 <- rowMeans(mtrx[, grps[[1]]]); 
  m2 <- rowMeans(mtrx[, grps[[2]]])/w;
  l2 <- log2(pmax(m2, min(m2[m2>0])/2)) - log2(pmax(m1, min(m1[m1>0])/2)); 
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps);
}