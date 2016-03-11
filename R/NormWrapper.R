######################################################################################
# Names of normalization methods
NormMethods<-function () c('NormLoess', 'NormAffyLoess');

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
  
  # return a list
  list(
    data=mtrx,
    method=mthd,
    parameters=args,
    results=do.call(mthd, all.args) # call the selected method with an argument list
  )
}

NormLoess<-function(mtrx, ref=c('mean', 'median', 'first', 'last'), ...) {
  
  ref<-tolower(ref)[1];
  
  if (ref=='median') x<-apply(mtrx, 1, median) else
    if (ref=='first') x<-mtrx[, 1] else 
      if (ref=='last') x<-mtrx[, ncol(mtrx)] else 
        x<-rowMeans(mtrx);
  
  apply(mtrx, 2, function(y) {
    z<-residuals(loess(y~x));
    z+x;
  }); 
}


NormAffyLoess<-function(mtrx, log.it=FALSE, ...) {
  library(affy);
  normalize.loess(mtrx, log.it=log.it);
}

NormAffyQspline<-function(mtrx, ...) {
  library(affy);
  d<-normalize.qspline(mtrx);
  dimnames(d)<-dimnames(mtrx);
  d;
}

NormAffyConstant<-function(mtrx, ref=c('mean', 'median', 'first', 'last'), ...) {
  
  ref<-tolower(ref)[1];
  
  if (ref=='median') x<-apply(mtrx, 1, median) else
    if (ref=='first') x<-mtrx[, 1] else 
      if (ref=='last') x<-mtrx[, ncol(mtrx)] else 
        x<-rowMeans(mtrx);
  
  normalize.constant(mtrx, x);
}


