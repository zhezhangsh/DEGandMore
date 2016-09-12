# BADER
# https://www.bioconductor.org/packages/3.3/bioc/html/BADER.html
# BADER {BADER}
DeBader <- function(mtrx, grps, paired=FALSE, normalized=FALSE, reps = 10000, cl = 4, ...) {
  require(DEGandMore);
  require(BADER);
  require(parallel);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by BADER; performing unpaired test instead.\n");
  
  design <- factor(rep(names(grps), sapply(grps, length)));  
  
  cl <- max(1, cl); 
  rp <- ceiling(reps/cl);
  res <- mclapply(1:cl, function(i) {
    BADER::BADER(mtrx, design, sizeFactors = TRUE, reps = rp, burn = 1000, printEvery = 10000, saveEvery = 10000, mode = 'minimal')
  }, mc.cores = cl);

  m1 <- exp(rowMeans(sapply(res, function(res) res$logMeanA)));
  m2 <- exp(rowMeans(sapply(res, function(res) res$logMeanB)));
  l2 <- log2(exp(rowMeans(sapply(res, function(res) res$logFoldChange)))); 
  pv <- 1 - 2 * (rowMeans(sapply(res, function(res) res$diffProb)) - 0.5); 
  pv[pv==0 & !is.na(pv)] <- 1/reps;
  qv <- p.adjust(pv, method='BH');

  m1[rowSums(mtrx[, grps[[1]]])==0] <- 0;
  m2[rowSums(mtrx[, grps[[2]]])==0] <- 0;
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps)[2:1], collapse='-'), 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, bader=res);
}

