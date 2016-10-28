# Clean up the statistic table of DE analysis
OutputDe <- function(s, mtrx, grps, pvalue.offset.ratio=0) {
  colnames(s)[1:6] <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  s[is.na(s[,1]), 1] <- 0; 
  s[is.na(s[,2]), 2] <- 0; 
  s[is.na(s[,4]), 4] <- 0; 
  s[is.na(s[,5]), 5] <- 1; 
  s[is.na(s[,6]), 6] <- 1; 
  s[,3] <- s[,2]-s[,1]; 
  s[,5] <- pmin(1, s[, 5]);
  s[,6] <- pmin(1, s[, 6]);
  
  if (pvalue.offset.ratio > 0) {
    s[s[, 5]==0, 5] <- pvalue.offset.ratio*min(s[s[,5]>0, 5]); 
  }
  
  s; 
}