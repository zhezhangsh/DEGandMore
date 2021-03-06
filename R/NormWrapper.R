######################################################################################
##################################################################################################################
# Names of normalization methods
NormMethods<-function () {
  c("NormAffyConstant",    # normalize.constant() function of the affy package, using constant scaling factors
    "NormAffyLoess",       # normalize.loess() function of the affy package, fitting loess normalization
    "NormAffyQspline",     # normalize.qspline() function of the affy package, using quantiles to fit cubic splines
    "NormCyclicLoess",     # normalizeCyclicLoess() function of the limma package
    "NormDESeq",           # estimateSizeFactors() function of the DESeq2 package for RNA-seq, median of ratio
    "NormFPKM",            # fragment per kilobases per million reads, RNA-seq only
    "NormLinear",          # lm() function of the base package, rescale by fitting linear regression
    "NormLoess",           # rescale by fitting loess regression to a reference sample
    "NormMedian",          # rescale by median of non-zero genes
    "NormQQ",              # quantile-quantile normalization
    "NormRLE",             # calcNormFactors() function of the edgeR package, using the relative log (RLE) method
    "NormRlog",            # rlog() function of the DESeq2 function
    "NormTMM",             # calcNormFactors() function of the edgeR package, using the weighted trimmed mean of M-values (TMM) method
    "NormTotalCount",      # rescale by total read count, RNA-seq only
    "NormTPM",             # transcripts per million, RNA-seq only
    "NormUpperQuantile",   # rescale by upper quantile of non-zero genes
    "NormVST")             # varianceStabilizingTransformation() function of the DESeq2 function
}
##################################################################################################################

NormWrapper <- function(mtrx, mthd=NormMethods[1], args=list()) {
  # mtrx  A numeric matrix of gene expression data. Rows are samples
  # mthd  Name of the method to use for normalization
  # args  Name:value pairs of specific arguments to selected method 
  
  library(DEGandMore);
  
  if (!(mthd %in% NormMethods())) stop("No method named: ", mthd, "\n");
  mtrx<-as.matrix(mtrx);
  
  # Full argument list
  all.args<-append(list(mtrx=mtrx), args);
  
  do.call(mthd, all.args) # call the selected method with an argument list
}

####################################################################################
# Rescale by fitting linear regression
NormLinear<-function(mtrx, ref=c('mean', 'median', 'first', 'last'), trim=c(0, 1)) {
  # ref:    How the reference X will be selected
  # trim:   Lower and upper fraction to trim the reference X as the training data for linear model
  ref <- tolower(ref)[1];
  if (ref[1]=='median') x <- apply(mtrx, 1, median) else
    if (ref[1]=='first') x <- mtrx[, 1] else 
      if (ref[1]=='last') x <- mtrx[, ncol(mtrx)] else 
        x <- rowMeans(mtrx);
      
  trim <- sort(trim); 
  if (trim[1]>0 | trim[2]<1) {
    y   <- sort(x); 
    lo  <- y[max(1, floor(trim[1]*nrow(mtrx)))];
    hi  <- y[min(nrow(mtrx), ceiling(trim[2]*nrow(mtrx)))]; 
    ind <- which(x>=lo & x<=hi);
  } else ind <- 1:nrow(mtrx);
  
  d <- apply(mtrx, 2, function(y) {
    mdl <- lm(y[ind] ~ x[ind]); 
    cof <- coefficients(mdl);
    prd <- (x - cof[1]) / cof[2];
    y - prd + x;
  }); 
  
  dimnames(d) <- dimnames(mtrx); 
  
  d; 
}

####################################################################################
# Rescale by fitting loess regression
NormLoess<-function(mtrx, ref=c('mean', 'median', 'first', 'last'), thread=4, degree=1, seed=1:nrow(mtrx), ...) {
  # degree:  the degree of the polynomials to be used, normally 1 or 2.
  
  ref <- tolower(ref)[1];
  if (ref[1]=='median') x <- apply(mtrx, 1, median) else
    if (ref[1]=='first') x <- mtrx[, 1] else 
      if (ref[1]=='last') x <- mtrx[, ncol(mtrx)] else 
        x <- rowMeans(mtrx, na.rm=TRUE);

  rnk  <- rank(x, ties.method = 'random'); 
  d0   <- lapply(1:ncol(mtrx), function(i) mtrx[, i]); 
  s0   <- c(which(rnk<=1000), which(rnk>=(nrow(mtrx)-999)));
  seed <- union(seed, s0); 
  sp   <- min(0.9, round(2000/length(seed)+0.05, 1)); 
      
  d <- parallel::mclapply(d0, function(y) {
    dat <- data.frame(x=x, y=y); 
    mdl <- loess(y ~ x, degree = degree, span = sp, data=dat[seed, , drop=FALSE]);
    prd <- predict(mdl, newdata=data.frame(x=x)); 
    res <- y - prd;
    res + x;
  }, mc.cores = thread);
  d <- do.call('cbind', d);
  dimnames(d) <- dimnames(mtrx);
  
  d;
}

####################################################################################
NormTotalCount <- function(mtrx, ref=c('mean', 'median', 'first', 'last')) {
  ref<-tolower(ref)[1];
  
  if (ref[1]=='median') x<-apply(mtrx, 1, median) else
    if (ref[1]=='first') x<-mtrx[, 1] else 
      if (ref[1]=='last') x<-mtrx[, ncol(mtrx)] else 
        x<-rowMeans(mtrx, na.rm=TRUE);
      
  x[is.na(x)] <- rowMeans(mtrx, na.rm=TRUE)[is.na(x)]; 
  
  apply(mtrx, 2, function(c) c/(mean(c)/mean(x, na.rm=TRUE)));
}

####################################################################################
NormQQ <- function(mtrx, ref=c('mean', 'median', 'first', 'last')) {
  ref<-tolower(ref)[1];
  
  if (ref[1]=='median') x<-apply(mtrx, 1, median) else
    if (ref[1]=='first') x<-mtrx[, 1] else 
      if (ref[1]=='last') x<-mtrx[, ncol(mtrx)] else 
        x<-rowMeans(mtrx, na.rm=TRUE);
  x <- rev(sort(x));
  
  d <- apply(mtrx, 2, function(y) {
    names(y) <- 1:length(y); 
    z <- y[y > x[length(x)]]; 
    z[1:length(z)] <- x[length(z)-rank(z, ties.method='random')+1]; 
    y[names(z)] <- z; 
    as.vector(y); 
  });  
  
  rownames(d) <- rownames(mtrx); 
  d;
}

####################################################################################
NormUpperQuantile <- function(mtrx, ref=c('mean', 'median', 'first', 'last')) {
  ref<-tolower(ref)[1];
  
  if (ref[1]=='median') x<-apply(mtrx, 1, median) else
    if (ref[1]=='first') x<-mtrx[, 1] else 
      if (ref[1]=='last') x<-mtrx[, ncol(mtrx)] else 
        x<-rowMeans(mtrx, na.rm=TRUE);
  
  x[is.na(x)] <- rowMeans(mtrx, na.rm=TRUE)[is.na(x)]; 
  y <- x[x>0]; 
  y <- y[y<=quantile(y)[4]]; 
  z <- which(rownames(mtrx) %in% names(y)); 

  apply(mtrx, 2, function(c) c/(mean(c[z])/mean(x[z], na.rm=TRUE)));
}

####################################################################################
NormMedian <- function(mtrx, ref=c('mean', 'median', 'first', 'last')) {
  ref<-tolower(ref)[1];
  
  if (ref[1]=='median') x<-apply(mtrx, 1, median) else
    if (ref[1]=='first') x<-mtrx[, 1] else 
      if (ref[1]=='last') x<-mtrx[, ncol(mtrx)] else 
        x<-rowMeans(mtrx, na.rm=TRUE);
      
  x[is.na(x)] <- rowMeans(mtrx, na.rm=TRUE)[is.na(x)]; 
  y <- median(x[x>0], na.rm=TRUE); 

  apply(mtrx, 2, function(c) c/(median(c[!is.na(x) & x>0])/y));
}

####################################################################################
# FPKM
NormFPKM <- function(mtrx, len) {
  len[is.na(len)] <- median(len); 
  apply(mtrx, 2, function(x) x/sum(x)/len*10^9); 
}

####################################################################################
# TPM
NormTPM <- function(mtrx, len) {
  len[is.na(len)] <- median(len); 
  
  d <- apply(mtrx, 2, function(x) x/len*1000); 
  apply(d, 2, function(d) d/(sum(d)/10^6)); 
}

####################################################################################
# DESeq2
NormDESeq <- function(mtrx) {
  
  require(DESeq2); 
  
  suppressMessages({
    dds <- DESeqDataSetFromMatrix(mtrx, colData = DataFrame(factor(rep(1, ncol(mtrx)))), design = ~1);
    dds <- estimateSizeFactors(dds);
    dds <- estimateDispersions(dds);
  }); 

  DESeq2::counts(dds, normalized = TRUE); 
}

####################################################################################
# Regularized log transformation in DESeq2
NormRlog <- function(mtrx) {
  require(DESeq2); 
  
  d <- rlog(mtrx);
  dimnames(d) <- dimnames(mtrx); 
  
  d; 
}

####################################################################################
# Variance stabilizing transformation in DESeq2
NormVST <- function(mtrx) {
  require(DESeq2); 
  varianceStabilizingTransformation(mtrx); 
}

####################################################################################
# TMM in edgeR
NormTMM <- function(mtrx, ref=c('mean', 'median', 'first', 'last')) {
  require(edgeR);
  
  ref<-tolower(ref)[1];
  
  if (ref[1]=='median') x<-apply(mtrx, 1, median) else
    if (ref[1]=='first') x<-mtrx[, 1] else 
      if (ref[1]=='last') x<-mtrx[, ncol(mtrx)] else 
        x<-rowMeans(mtrx, na.rm=TRUE);
  
  dge <- DGEList(cbind(x, mtrx)); 
  dge <- calcNormFactors(dge, refColumn = 1, method = 'TMM'); 
  dge <- estimateDisp(dge);
  dge <- estimateCommonDisp(dge);
  dge <- estimateTagwiseDisp(dge);
  
  d <- dge@.Data[[11]][, -1];
  d;
}

####################################################################################
# RLE in edgeR
NormRLE <- function(mtrx, ref=c('mean', 'median', 'first', 'last')) {
  require(edgeR);
  
  ref<-tolower(ref)[1];
  
  if (ref[1]=='median') x<-apply(mtrx, 1, median) else
    if (ref[1]=='first') x<-mtrx[, 1] else 
      if (ref[1]=='last') x<-mtrx[, ncol(mtrx)] else 
        x<-rowMeans(mtrx, na.rm=TRUE);
      
  dge <- DGEList(cbind(x, mtrx)); 
  dge <- calcNormFactors(dge, refColumn = 1, method = 'RLE'); 
  dge <- estimateCommonDisp(dge);
  dge <- estimateTagwiseDisp(dge);
      
  d <- dge@.Data[[4]][, -1];
  d;
}

####################################################################################
# normalize.loess {affy}
# http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/affy/html/normalize.loess.html
NormAffyLoess<-function(mtrx, log.it=FALSE, ...) {
  require(affy);
  normalize.loess(mtrx, log.it=log.it);
}

####################################################################################
# normalize.qspline {affy}
# http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/affy/html/normalize.qspline.html
NormAffyQspline<-function(mtrx, ...) {
  require(affy);
  d<-normalize.qspline(mtrx);
  dimnames(d)<-dimnames(mtrx);
  d;
}

####################################################################################
# normalize.constant {affy}
# http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/affy/html/normalize.constant.html
NormAffyConstant<-function(mtrx, ref=c('mean', 'median', 'first', 'last'), ...) {
  require(affy);
  
  ref<-tolower(ref)[1];
  if (ref[1]=='median') x<-apply(mtrx, 1, median) else
    if (ref[1]=='first') x<-mtrx[, 1] else 
      if (ref[1]=='last') x<-mtrx[, ncol(mtrx)] else 
        x<-rowMeans(mtrx);
  
  normalize.constant(mtrx, x);
}

####################################################################################
# normalizeCyclicLoess {limma}
# http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/normalizeCyclicLoess.html
NormCyclicLoess <- function(mtrx) {
  require(limma);
  
  normalizeCyclicLoess(mtrx, method='pairs')
}
####################################################################################
