# Empirical Bayesian analysis of patterns of differential expression in count data
# baySeq package: http://bioconductor.org/packages/release/bioc/html/baySeq.html
DeBaySeq <- function(mtrx, grps, paired=FALSE, samplesize=10000, cl=2) {
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
  
  n   <- sapply(grps, length); 
  ncl <- cl[1];
  cl  <- NULL; 
  if (is.numeric(ncl[1])) {
    try(require(snow));
    try(cl <- makeCluster(max(1, round(ncl)), "SOCK")); 
  }
  
  if (paired & n[1]==n[2]) {
    cd <- new("countData", data = array(c(mtrx[, grps[[1]]], mtrx[, grps[[2]]]), dim = c(nrow(mtrx), n[1], 2)),
              replicates = c(1:n[1]), groups = list(NDE = rep(1, n[1]), DE = 1:n[1]), densityFunction = bbDensity); 
    libsizes(cd) <- getLibsizes(cd); 
    cd  <- getPriors(cd, samplesize = samplesize, cl = cl);
    cd  <- getLikelihoods(cd, pET = 'BIC', nullData = TRUE, cl = cl); 
    res <- topCounts(cd, 1, number = nrow(mtrx), normaliseData = TRUE); 
    pv  <- 1 - res$Likelihood; 
    nom <- lapply(1:n[1], function(x) {
      x <- as.vector(res[, x+1]); 
      x <- strsplit(x, ':'); 
      x <- do.call('rbind', x); 
      apply(x, 2, as.numeric);
    });
    m1  <- rowMeans(sapply(nom, function(x) x[, 1])); 
    m2  <- rowMeans(sapply(nom, function(x) x[, 2])); 
    rnm <- rownames(mtrx)[res[, 1]]; 
  } else {
    rp  <- rep(names(grps), n); 
    nde <- rep(1, ncol(mtrx)); 
    de  <- rep(1:2, n); 
    gp  <- list(NDE = nde, DE = de);
    cd  <- new("countData", data = mtrx, replicates = rp,  groups = gp); 
    libsizes(cd) <- getLibsizes(cd); 
    cd  <- getPriors.NB(cd, samplesize=samplesize, estimation = "QL", cl=cl);
    cd  <- getLikelihoods(cd, cl = cl, bootStraps = 1, verbose = FALSE); 
    res <- topCounts(cd, 2, number=nrow(mtrx), normaliseData = TRUE); 
    pv  <- 1- res$Likelihood; 
    m1  <- rowMeans(res[, colnames(mtrx)[grps[[1]]]]); 
    m2  <- rowMeans(res[, colnames(mtrx)[grps[[2]]]]); 
    rnm <- rownames(mtrx)[res[, 1]]; 
  }; 
  
  if (!is.null(cl)) stopCluster(cl); 

  qv  <- p.adjust(pv, method='BH');
  l2 <- CalculateCountLog2FC(m2, m1, mtrx, grps);
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rnm; 
  
  list(stat=s[rownames(mtrx), ], group=grps, cd=cd);
}