DeBBSeq <- function(mtrx, grps, paired=FALSE) {
  require(DEGandMore);
  require(BBSeq); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  n <- sapply(grps, length);

  if (paired) warning("Paired test of BBSeq is disabled due to bug in source code.\n");
  
#   if (paired & n[1]==n[2]) {
#     a1 <- rep(0:1, each=n[1]); 
#     a2 <- rep(1:n[1], 2);
#     X <- X.builder(cbind(a1, a2), c(1, 0));
#     
#     libsize <- lib.size(mtrx); 
#     like = like.matrix(mtrx, libsize, c(1, 0), cbind(a1, a2)); 
#   }
  
  X <- cbind(1, rep(0:1, n)); 
  output <- free.estimate(mtrx, X); 
  out.model <- constrained.estimate(mtrx, X, gn=3, output$betahat.free, psi.free=output$psi.free);
  flag <- outlier.flag(mtrx); 

  pv <- out.model$p.model; 
  qv <- p.adjust(pv, method='BH');
  
  l2 <- log2(exp(out.model$betahat.model[, 2])); 
  tt <- sum(rowMeans(mtrx));
  m0 <- out.model$betahat.model[, 1];
  m1 <- m2 <- exp(m0 + log(tt));
  fc <- 2^(l2/2); 
  m1 <- m1/fc;
  m2 <- m2*fc;
  
  m1[is.na(pv)] <- 0;
  m2[is.na(pv)] <- 0;
  l2[is.na(pv)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv, flag);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR', 'Flag');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps);
}