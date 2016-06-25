
######################################################################################
# Names of DE methods
DeMethods <- function () c('DeT', 'DeSam', 'DeRankP', 'DeDeSeq', 'DeEdgeR', 'DeVoomLimma');

# DeT         Student's T test
# DeSam       SAM (significance analysis of microarray) test
# DeRankP     Rank product test
# DeDeSeq     DeSeq test for RNA-seq data
# DeEdgeR     EdgeR test for RNA-seq data
# DeVoomLimma Voom normalization followed by Limma for DE
######################################################################################

# Whether the method is applicable to count data only
DeCountMethods <- function () c('DeT'=FALSE, 'DeSam'=FALSE, 'DeRankP'=FALSE, 'DeDeSeq'=TRUE, 'DeEdgeR'=TRUE, 'DeVoomLimma'=TRUE);

DeWrapper <- function(mtrx, grps, mthd=DeMethods()[1], paired=FALSE, logged=TRUE, args=list()) {
    # mtrx    A numeric matrix of gene expression data. Rows are unique genes and columns include 2 groups of samples to be compared
    # grps    A list of two vectors, each vector has the column indexes (numeric vectors) or column names (character vectors) of a group; vectors are named by group names  
    # mthd    Name of the method to use for differential expression analysis
    # paired  Whether it's a paired test; if TRUE the two groups must have the same number of samples and the samples must be ordered to match each other
    # logged  Whether data is log2-transformed. It will affect how the LogFC values, but not the p values, will be calculated. It won't affect methods using RNA-seq read count data either.
    # args    Name:value pairs of specific arguments to selected method 
  
    require(DEGandMore);
  
    if (!(mthd %in% DeMethods())) stop("DE method '", mthd, "' not available.\n") else mthd <- mthd[1];
    
    # standardize data types of inputs to all methods
    mtrx <- as.matrix(mtrx);
    grps <- lapply(grps[1:2], function(g) 
      if (class(g) == 'character') which(colnames(mtrx) %in% g) else g[g>0 & g<=ncol(mtrx)]);
    
    # Give error if no qualified samples for comparison
    if (min(sapply(grps, length)) == 0) 
      stop("No samples found in data matrix for comparison.\n");
    if (paired & length(grps[[1]])!=length(grps[[2]])) 
      stop("Numbers of samples in groups are unequal for paired test.\n");
    
    grp0 <- colnames(mtrx)[grps[[1]]];
    grp1 <- colnames(mtrx)[grps[[2]]];
    
    # Full argument list
    all.args <- list(mtrx=mtrx, grps=grps, paired=paired);
    if (DeCountMethods()[mthd]) logged <- FALSE else all.args <- append(all.args, logged); 
    all.args <- append(all.args, args);
    
    res <- do.call(mthd, all.args); # call the selected method with an argument list
    res$stat <- res$stat[rownames(mtrx), , drop=FALSE];
    
    # return a list
    list(
      data       = mtrx,
      group0     = grp0,
      group1     = grp1,
      method     = mthd,
      paired     = paired,
      logged     = logged,
      parameters = args,
      results    = res  
    )
  }
