
######################################################################################
# Names of DE methods
DeMethods <- function () c('DeT', 'DeSam', 'DeRankP', 'DeDeSeq2', 'DeEdgeR', 'DeVoomLimma');

# DeT         Student's T test
# DeSam       SAM (significance analysis of microarray) test
# DeRankP     Rank product test
# DeDeSeq2    DeSeq test for RNA-seq data
# DeEdgeR     EdgeR test for RNA-seq data
# DeVoomLimma Voom normalization followed by Limma for DE
######################################################################################

# Whether the method is applicable to count data only
DeCountMethods <- function () c('DeT'=FALSE, 'DeSam'=FALSE, 'DeRankP'=FALSE, 'DeDeSeq2'=TRUE, 'DeEdgeR'=TRUE, 'DeVoomLimma'=TRUE);

DeWrapper <- function(mtrx, grps, mthd=DeMethods()[1], paired=FALSE, logged=FALSE, args=list()) {
    # mtrx    A numeric matrix of gene expression data. Rows are unique genes and columns include 2 groups of samples to be compared
    # grps    A list of two vectors, each vector has the column indexes (numeric vectors) or column names (character vectors) of a group; vectors are named by group names  
    # mthd    Name of the method to use for differential expression analysis
    # paired  Whether it's a paired test; if TRUE the two groups must have the same number of samples and the samples must be ordered to match each other
    # args    Name:value pairs of specific arguments to selected method 
  
    library(DEGandMore);
  
    if (!(mthd %in% DeMethods())) stop("DE method '", mthd, "' not available.\n") else mthd <- mthd[1];
    
    # standardize data types of inputs to all methods
    mtrx <- as.matrix(mtrx);
    grps <- lapply(grps[1:2], function(g) 
      if (class(g) == 'character') which(colnames(mtrx) %in% g) else g[g>0 & g<=ncol(mtrx)]);
    
    # Give error if no qualified samples for comparison
    if (min(sapply(grps, length)) == 0) 
      stop("No samples found in data matrix for comparison.\n");
    if (pair & length(grps[[1]])!=length(grps[[2]])) 
      stop("Numbers of samples in groups are unequal for paired test.\n");
    
    grp0 <- colnames(mtrx)[grps[[1]]];
    grp1 <- colnames(mtrx)[grps[[2]]];
    
    # Full argument list
    all.args <- list(mtrx=mtrx, grps=grps, paired=paired);
    if (DeCountMethods()[mthd[1]])
    all.args <- append(all.args, args);
    
    # return a list
    list(
      data       = mtrx,
      group0     = grp0,
      group1     = grp1,
      method     = mthd,
      paired     = pair,
      parameters = args,
      results    = do.call(mthd, all.args) # call the selected method with an argument list
    )
  }
