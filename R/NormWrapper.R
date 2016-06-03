######################################################################################
# Names of normalization methods
NormMethods<-function () 
  c("NormAffyConstant", "NormAffyLoess", "NormAffyQspline", "NormDESeq", "NormFPKM", "NormLoess", "NormMedian", "NormRLE", "NormTMM", 
    "NormTotalCount", "NormTPM", "NormUpperQuantile", "NormWrapper", 'NormLoess', 'NormAffyLoess', 'NormAffyQspline');

# NormLoess           Self-implemented Loess normalization
# NormAffyLoess       normalize.loess function of Affy package
# NormAffyQspline     normalize.qsline function of Affy package
# NormAffyConstant    normalize.constant function of Affy package 
######################################################################################

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
# DESeq
NormDESeq <- function(mtrx) {
  require(DESeq); 
  
  cds <- newCountDataSet(mtrx, rep('A', ncol(mtrx)));
  cds <- estimateSizeFactors(cds); 
  f <- sizeFactors(cds);
  
  d <- sapply(1:ncol(mtrx), function(i) mtrx[, i]/f[i]); 
  colnames(d) <- colnames(mtrx); 
  
  d;
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
  
  f <- calcNormFactors(cbind(x, mtrx), refColumn = 1, method = 'RLE')[-1];
  d <- sapply(1:ncol(mtrx), function(i) mtrx[, i]/f[i]); 
  colnames(d) <- colnames(mtrx); 
  
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
      
  f <- calcNormFactors(cbind(x, mtrx), refColumn = 1, method = 'RLE')[-1];
  d <- sapply(1:ncol(mtrx), function(i) mtrx[, i]/f[i]); 
  colnames(d) <- colnames(mtrx); 
      
  d;
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
NormLoess<-function(mtrx, ref=c('mean', 'median', 'first', 'last'), thread=4, ...) {
  
  ref<-tolower(ref)[1];
  if (ref[1]=='median') x<-apply(mtrx, 1, median) else
    if (ref[1]=='first') x<-mtrx[, 1] else 
      if (ref[1]=='last') x<-mtrx[, ncol(mtrx)] else 
        x<-rowMeans(mtrx);
  
  d0<-lapply(1:ncol(mtrx), function(i) mtrx[, i]); 
  
  d<-parallel::mclapply(d0, function(y) {
    z<-residuals(loess(y~x));
    z+x;
  }, mc.cores = thread); 

  d<-do.call('cbind', d);
  colnames(d)<-colnames(mtrx); 
  
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
