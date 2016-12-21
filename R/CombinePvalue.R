# Combine p values
CombinePvalueMethods <- function() {
  c(Fisher='fisher', Simes='simes', Bonferroni='bonferroni', Maximum='max', Minimum='min', Average='average'); 
}

CombinePvalue <- function(pv, mthd=c('fisher', 'simes', 'bonferroni', 'max', 'min', 'average'), 
                          offset.ratio=0.5, normalize=TRUE, adjust.fisher=TRUE) {
  # pv      A matrix of p values to combine the columns
  # mthd    Method to combined p values
            # 'fisher': Fisher's meta analysis, default method
            # 'simes': Simes' method
            # 'bonferroni': bonferroni adjustment by number of p values per gene
            # 'max': maximum of all p values
            # 'min': minimum of all p values
            # 'average': geometric mean
  # offset.ratio  How to replace p value of 0, the ratio to the minimal non-zero p value
  # normalize     If TRUE, QQ normalize the p value sets first 
  
  # pv[is.na(pv)] <- 1;
  pv[pv < 0] <- 0;
  pv[pv > 1] <- 1;
  mn <- min(pv[pv>0]); 
  pv[pv==0] <- mn*offset.ratio;
  pv[pv==0] <- mn; # minimum possible value
  
  if (normalize) {
    st <- apply(pv, 2, sort); 
    mn <- 10^rowMeans(log10(st)); 
    pv <- apply(pv, 2, function(p) mn[round(rank(p))]); 
  }
  
  mthd <- tolower(mthd[1]); 
  
  if (mthd == 'bonferroni') {
    p <- apply(pv, 1, function(p) length(n[!is.na(n)])*min(p, na.rm=TRUE));
    p <- pmin(1, p); 
  } else if (mthd == 'max') {
    p <- apply(pv, 1, function(p) max(p, na.rm=TRUE)); 
  } else if (mthd == 'min') {
    p <- apply(pv, 1, function(p) min(p, na.rm=TRUE)); 
  } else if (mthd == 'average') {
    p <- 10^rowMeans(log10(pv), na.rm=TRUE); 
  } else if (mthd == 'simes') {
    n <- ncol(pv); 
    p <- apply(pv, 1, function(p) {
      n <- length(p[!is.na(p)]); 
      s <- min(n*(sort(p)/(1:n)), na.rm=TRUE); 
    }); 
    p <- pmin(1, p); 
  } else {
    p <- apply(pv, 1, function(p) {
      n <- length(p[!is.na(p)]); 
      c <- -2*sum(log(p), na.rm=TRUE);
      p <- pchisq(c, 2*n, lower.tail = FALSE, log.p = TRUE); 
      if (adjust.fisher) p <- exp(p + log(n)) else p <- exp(p); 
      p; 
    }); 
    p <- pmin(1, p); 
  }
  
  p;  
}