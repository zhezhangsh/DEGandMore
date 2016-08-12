# Internal function that prepares data for differential expression analysis
PrepareDe <- function(mtrx, grps, pair) {
  
  # Groups
  grps <- lapply(grps[1:2], function(g) 
    if (class(g) == 'character') which(colnames(mtrx) %in% g) else g[g>0 & g<=ncol(mtrx)]);
  if (min(sapply(grps, length)) <= 1) stop("Not enough observation, require at least 2 samples in each group.\n");
  if (pair) if (length(grps[[1]]) != length(grps[[2]])) {
    warning('Warning: number of samples not equal between 2 groups; running unpaired test instead.\n'); 
    pair <- FALSE;
  }; 
  if (is.null(names(grps)[1])) names(grps)[1]<-'Group0';
  if (is.null(names(grps)[2])) names(grps)[2]<-'Group1';
  
  if (class(mtrx) != 'matrix') mtrx<-as.matrix(mtrx);
  e1 <- mtrx[, grps[[1]], drop=FALSE];
  e2 <- mtrx[, grps[[2]], drop=FALSE];
  d  <- cbind(e1, e2);
  
  grps[[1]] <- 1:ncol(e1);
  grps[[2]] <- (ncol(e1)+1):(ncol(e1)+ncol(e2)); 
  
  list(d, grps, pair);
}