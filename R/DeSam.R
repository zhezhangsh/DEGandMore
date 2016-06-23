# Perform a 2 group comparison with SAM method
DeSam<-function(mtrx, grps, paired=FALSE, logged=TRUE, nperm=100, ...) {
  # mtrx          A numeric matrix of gene expression data. Rows are unique genes and columns include 2 groups of samples to be compared
  # grps          A 2-vector list, each vector has the column indexes (numeric vectors) or column names (character vectors) of a group; vectors are named by group names  
  # paired  		  Paired test
  # logged        If the data has been log-transformed
  # nperm				  SAM parameter
  
  require(DEGandMore);
  require(samr);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  e1 <- mtrx[, grps[[1]], drop=FALSE];
  e2 <- mtrx[, grps[[2]], drop=FALSE];
  if (paired) {
    y <- rep(c(-1, 1), sapply(grps, length));
    type <- "Two class paired";
  } else {
    y <- rep(c(1, 2), sapply(grps, length));
    type <- "Two class unpaired";
  }
  
  data <- list(x=cbind(e1, e2), y=y, geneid=rownames(mtrx), genenames=rownames(mtrx), logged2=logged); 
  
  capture.output(sam<-samr(data, resp.type=type, nperms=nperm, ...))->x;
  p <- samr.pvalues.from.perms(sam$tt, sam$ttstar);
  capture.output(delta<-samr.compute.delta.table(sam))->x;
  x <- samr.compute.siggenes.table(sam, 0, data, delta, all.genes=T); 
  x <- rbind(x[[1]], x[[2]]);  
  rownames(x) <- as.vector(x[,2]);  
  x  <- x[rownames(mtrx), ]; 
  fc <- as.numeric(x[,7]); 
  q  <- as.numeric(x[,8])/100; 
  m1 <- rowMeans(e1); 
  m2 <- rowMeans(e2);
  nm <- paste(names(grps), collapse='-');
  if (logged) lgfc<-m2-m1 else lgfc<-log2(fc);
  lgfc[is.na(lgfc)]<-0; 
  s <- cbind(m1, m2, m2-m1, lgfc, pmin(1, pmax(0, p)), pmin(1, pmax(0, q)));
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps), collapse='-'), 'LogFC', 'Pvalue', 'FDR');
  
  list(stat=s[rownames(mtrx), ], group=grps, sam=sam, sig=x)
}