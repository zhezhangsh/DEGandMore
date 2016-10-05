# Empirical Bayesian analysis of patterns of differential expression in count data
# baySeq package: http://bioconductor.org/packages/release/bioc/html/baySeq.html
DeBaySeq <- function(mtrx, grps, paired=FALSE, normalization=c('None', 'DESeq', 'TMM', 'RLE', 'Median', 'UQ', 'TC', 'QQ'), 
                     samplesize=1000, bootStraps=2, cl=4) {
  # normalization methods:
    # TMM: Trimmed mean method of calcNormFactors() function of the edgeR package
    # RLS: Relative log method of calcNormFactors() function of the edgeR package
    # DESeq: estimateSizeFactors() function of the DESeq2 package for RNA-seq, median of ratio
    # Median: rescale by median of non-zero genes
    # UQ: rescale by upper quantile of non-zero genes
    # TC: rescale by total read count
    # QQ: quantile-quantile normalization
  
  require(DEGandMore);
  require(baySeq);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  norm <- tolower(normalization)[1]; 
  if (norm=='deseq')  mtrx <- NormDESeq(mtrx);
  if (norm=='tmm')    mtrx <- NormTMM(mtrx);
  if (norm=='rle')    mtrx <- NormRLE(mtrx);
  if (norm=='median') mtrx <- NormMedian(mtrx);
  if (norm=='uq')     mtrx <- NormUpperQuantile(mtrx);
  if (norm=='tc')     mtrx <- NormTotalCount(mtrx);
  if (norm=='qq')     mtrx <- NormQQ(mtrx);
  mtrx0 <- mtrx;
  mtrx  <- round(mtrx); 
  mtrx  <- mtrx[, c(grps[[1]], grps[[2]]), drop=FALSE]; 
  
  n <- sapply(grps, length); 
  
  ncl <- cl[1];
  cl  <- NULL; 
  if (is.numeric(ncl[1])) {
    try(require(snow));
    try(cl <- makeCluster(max(1, round(ncl)), "SOCK")); 
  }
  
  if (paired & n[1]==n[2]) {
    d  <- list(mtrx[, grps[[1]]], mtrx[, grps[[2]]]); 
    rp <- 1:n[1]; 
    gp <- list(NDE = rep(1, n[1])); 
    cd <- new("countData", data = d, replicates = rp,  groups = gp, densityFunction = bbDensity); 
    libsizes(cd) <- getLibsizes(cd); 
    cd  <- getPriors(cd, samplesize = samplesize, cl = cl);
    cd  <- getLikelihoods(cd, pET = 'BIC', nullData = TRUE, bootStraps = bootStraps, cl = cl); 
    p <- 1 - exp(cd@posteriors[, 1]); 
  } else {
    rp  <- rep(names(grps), n); 
    nde <- rep(1, ncol(mtrx)); 
    de  <- rep(1:2, n); 
    gp  <- list(NDE = nde, DE = de);
    cd  <- new("countData", data = mtrx, replicates = rp,  groups = gp); 
    libsizes(cd) <- getLibsizes(cd); 
    cd  <- getPriors.NB(cd, samplesize=samplesize, estimation = "QL", cl=cl);
    cd  <- getLikelihoods(cd, cl = cl, bootStraps = bootStraps, verbose = FALSE); 
    p <- 1 - exp(cd@posteriors[, 2]); 
  }; 
  
  if (!is.null(cl)) stopCluster(cl); 

  m1 <- rowMeans(mtrx0[, grps[[1]], drop=FALSE], na.rm=TRUE);
  m2 <- rowMeans(mtrx0[, grps[[2]], drop=FALSE], na.rm=TRUE);
  q  <- p.adjust(p, method='BH');
  fc <- CalculateCountLog2FC(m1, m2, mtrx, grps);
  
  s <- cbind(m1, m2, m2-m1, fc, p, q);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s[rownames(mtrx), ], group=grps, cd=cd);
}