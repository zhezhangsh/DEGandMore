# Paired/Unpaired Student's t test
# Equal variance
DeT <- function(mtrx, grps, paired=FALSE, logged=TRUE) {
  # mtrx    A numeric matrix of gene expression data. Rows are unique genes and columns include 2 groups of samples to be compared
  # grps    A 2-vector list, each vector has the column indexes (numeric vectors) or column names (character vectors) of a group; vectors are named by group names  
  # paired  Whether it's a paired test; if TRUE the two groups must have the same number of samples and the samples must be ordered to match each other
  # logged  Whether data is log2-transformed. It will affect how Log(FoldChange), but not the p values, will be calculated. It won't affect methods using RNA-seq read count data either.
  
    require(DEGandMore);

    prepared <- PrepareDe(mtrx, grps, paired);
    mtrx     <- prepared[[1]];
    grps     <- prepared[[2]];
    paired   <- prepared[[3]];
    
    # data subsets of the sample groups
    d1 <- mtrx[, grps[[1]], drop=FALSE];
    d2 <- mtrx[, grps[[2]], drop=FALSE];
    
    # sample sizes
    n1 <- apply(d1, 1, function(d) length(d[!is.na(d)]));
    n2 <- apply(d2, 1, function(d) length(d[!is.na(d)]));
    n  <- sapply(grps, length); 
    
    # group means
    m1 <- rowMeans(d1, na.rm=TRUE);
    m2 <- rowMeans(d2, na.rm=TRUE);
    diff <- m2-m1;
    
    if (logged) lgfc <- diff else {
      mn <- min(mtrx[mtrx>0]); 
      lgfc <- log2(pmax(mn, m2)) - log2(pmax(mn, m1)); 
    };
    
    # calculate variance while assuming unequal variance
    if (!paired | n[1]!=n[2]) {
      ss1 <- rowSums((d1-m1)^2, na.rm=TRUE)/(n1-1);
      ss2 <- rowSums((d2-m2)^2, na.rm=TRUE)/(n2-1);
      t   <- diff/sqrt(ss1/n1+ss2/n2); 
      df  <- n1+n2-2; # Degree of freedom
      t[is.na(t)] <- 0;
    } else {
      d   <- d2-d1; 
      n   <- apply(d, 1, function(d) length(d[!is.na(d)]));
      t   <- rowMeans(d, na.rm=TRUE)/(apply(d, 1, function(d) sd(d, na.rm=TRUE))/sqrt(n));
      df  <- n-1;
      t[is.na(t)] <- 0;
    }
     
    # calculated 2-sided p values
    p      <- rep(0.5, length(t));
    p[t>0] <- 1-pt(t[t>0], df=df[t>0]);
    p[t<0] <- pt(t[t<0], df=df[t<0]);
    p      <- p*2;
    
    # Return a statistic matrix
    stat <- cbind(m1, m2, diff, lgfc, p, p.adjust(p, method='BH'));
    rownames(stat) <- rownames(mtrx);
    colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');

    list(stat=stat[rownames(mtrx), ], group=grps);
  }