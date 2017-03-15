# Summary row statistics for clustering analysis: mean, range, and variance 
SummarizeStatCl <- function(d, group = colnames(d)) {
  # d       data matrix
  # group   the grouping of columns into sample gorups
  # nnn     number of nearest neighbor
  
  gp <- split(1:ncol(d), group)[as.character(unique(group))]; 
  ns <- sapply(gp, length);
  
  if (max(ns) == 1) { # Columns are not grouped
    m <- rowMeans(d, na.rm=TRUE);
    r <- apply(d, 1, function(d) range(d, na.rm=TRUE));
    v <- apply(d, 1, function(d) sd(d, na.rm=TRUE)); 
    s <- cbind(Mean=m, Min=r[1, ], Max=r[2, ], SD=v);
    o <- list(replicates=FALSE, stat=s); 
  } else {
    mm <- rowMeans(d, na.rm=TRUE);
    ms <- sapply(gp, function(g) rowMeans(d[, g, drop=FALSE]));
    rg <- apply(ms, 1, function(d) range(d, na.rm=TRUE));

    st <- rowSums((d-mm)^2); # Total sum of squares
    sw <- rowSums(sapply(gp, function(g) {
      d0 <- d[, g, drop=FALSE];
      rowSums((d0-rowMeans(d0))^2);
    })); 
    sb <- st-sw;
    fs <- (sb/(length(gp)-1))/(sw/(ncol(d)-length(gp)));
    s <- cbind(Mean=mm, Min=rg[1, ], Max=rg[2, ], Fvalue=fs);
    o <- list(replicates=TRUE, stat=s); 
  }; 

  o; 
}