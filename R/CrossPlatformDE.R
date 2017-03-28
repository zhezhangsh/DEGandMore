# Make the differential expression measured by different platform (ex. microarray vs. RNA-seq) more compatible to each other
# by adjusting mean and variance
CrossPlatformDE <- function(d0, d1, ref0=NA, ref1=NA) {
  # d0    The reference data matrix
  # d1    The data matrix to be adjusted
  # ref0  Choose a subset of columns in d0 as the reference; use all columns if NA
  # ref1  Choose a subset of columns in d1 for estimating parameters; use all columns if NA
  
  if (identical(NA, ref0)) n0 <- ncol(d0) else n0 <- length(ref0);
  if (identical(NA, ref1)) n1 <- ncol(d1) else n1 <- length(ref1);
  
  calV <- function(d) sqrt(rowSums(apply(d, 2, function(x) x-rowMeans(d))^2)/ncol(d));
  
  if (identical(NA, ref0)) {
    m0 <- rowMeans(d0);
    v0 <- calV(d0);
  } else {
    m0 <- rowMeans(d0[, ref0]);
    v0 <- calV(d0[, ref0]);
  }
  
  if (identical(NA, ref1)) m1 <- rowMeans(d1) else m1 <- rowMeans(d1[, ref1]);
  
  adj2 <- adj1 <- apply(d1, 2, function(x) x-(m1-m0)); 
  
  for (i in 1:nrow(adj1)) {
    v <- v0[i]; 
    if (identical(NA, ref1)) x <- adj1[i, ] else x <- adj1[i, ref1];
    y <- x-mean(x);
    z <- v/sqrt(sum(y^2)/length(y)); 
    adj2[i, ] <- z*(adj1[i, ]-mean(x)) + mean(x); 
  }
  adj2; 
}; 