DeT <- function(mtrx, grps, paired=FALSE, logged=TRUE) {
    # mtrx          A numeric matrix of gene expression data. Rows are unique genes and columns include 2 groups of samples to be compared
    # grps          A 2-vector list, each vector has the column indexes (numeric vectors) or column names (character vectors) of a group; vectors are named by group names  
    # 
    # column index in range
    
    library(DEGandMore);
  
    grps<-lapply(grps, function(g) if (class(g) == 'character') which(colnames(mtrx) %in% g) else g[g>0 & g<=ncol(mtrx)]);
    
    # use default group names if not given
    if (is.null(names(grps)[1])) names(grps)[1]<-'Group0';
    if (is.null(names(grps)[2])) names(grps)[2]<-'Group1';
    
    # Give error if no qualified samples for comparison
    if (min(sapply(grps, length)) <= 1) stop("No enough observation for a 2-sample t test.");
    
    if (class(mtrx) != 'matrix') mtrx<-as.matrix(mtrx);
    
    # data subsets of the sample groups
    d1<-mtrx[, grps[[1]], drop=FALSE];
    d2<-mtrx[, grps[[2]], drop=FALSE];
    
    # sample sizes
    n1<-ncol(d1);
    n2<-ncol(d2);
    
    # group means
    m1<-rowMeans(d1, na.rm=TRUE);
    m2<-rowMeans(d2, na.rm=TRUE);
    diff<-m2-m1;
    
    if (logged) lgfc<-diff else lgfc<-log2(m2/pmax(m1, min(m1[m1>0])));
    
    # calculate variance while assuming unequal variance
    if (!paired | n1!=n2) {
      ss1<-rowSums((d1-m1)^2, na.rm=TRUE)/(n1-1);
      ss2<-rowSums((d2-m2)^2, na.rm=TRUE)/(n2-1);
      s<-sqrt(ss1/n1+ss2/n2); # pooled variance
      t<-diff/s; # T statistic
      t[is.na(t)]<-0;
      df<-(ss1/n1+ss2/n2)^2/((ss1/n1)^2/(n1-1)+(ss2/n2)^2/(n2-1)); # Degree of freedom
    } else {
      d<-d2-d1; 
      n<-apply(d, 1, function(d) length(d[!is.na(d)]));
      t<-rowMeans(d, na.rm=TRUE)/(apply(d, 1, function(d) sd(d, na.rm=TRUE))/sqrt(n));
      t[is.na(t)]<-0;
      df<-n-1;
    }
     
    # calculated 2-sided p values
    p<-rep(0.5, length(t));
    p[t>0]<-1-pt(t[t>0], df=df[t>0]);
    p[t<0]<-pt(t[t<0], df=df[t<0]);
    p<-p*2;
    
    # Return a statistic matrix
    stat<-cbind(m1, m2, diff, lgfc, p, p.adjust(p, method='BH'));
    rownames(stat)<-rownames(mtrx);
    colnames(stat)<-c(paste('Mean', names(grps), sep='_'), paste(names(grps), collapse='-'), 'LogFC', 'Pvalue', 'FDR');
    list(stat=stat[rownames(mtrx), ], group=grps);
  }