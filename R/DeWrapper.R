
######################################################################################
# Names of DE methods
DeMethods<-function () c('DeT', 'DeSam', 'DeRankP', 'DeDeSeq2', 'DeEdgeR', 'DeVoomLimma');

# DeT         Student's T test
# DeSam       SAM (significance analysis of microarray) test
# DeRankP     Rank product test
# DeDeSeq2    DeSeq test for RNA-seq data
# DeEdgeR     EdgeR test for RNA-seq data
# DeVoomLimma Voom normalization followed by Limma for DE
######################################################################################

DeWrapper <- function(mtrx, grps, mthd=DeMethods()[1], args=list()) {
    # mtrx  A numeric matrix of gene expression data. Rows are unique genes and columns include 2 groups of samples to be compared
    # grps  A list of vectors, each vector has the column indexes (numeric vectors) or column names (character vectors) of a group; vectors are named by group names  
    # mthd  Name of the method to use for differential expression analysis
    # args  Name:value pairs of specific arguments to selected method 
  
    library(DEGandMore);
  
    if (!(mthd %in% DeMethods())) stop("No method named: ", mthd, "\n");
    
    # standardize data types of inputs to all methods
    mtrx<-as.matrix(mtrx);
    grps<-lapply(grps, function(g) if (class(g) == 'character') which(colnames(mtrx) %in% g) else g[g>0 & g<=ncol(mtrx)]);
    
    # Give error if no qualified samples for comparison
    if (min(sapply(grps, length)) == 0) stop("No samples for comparison");
    
    # Full argument list
    all.args<-list(mtrx=mtrx, grps=grps);
    all.args<-append(all.args, args);
    
    # return a list
    list(
      data=mtrx,
      group0=colnames(mtrx)[grps[[1]]],
      group1=colnames(mtrx)[grps[[2]]],
      method=mthd,
      parameters=args,
      results=do.call(mthd, all.args) # call the selected method with an argument list
    )
  }
