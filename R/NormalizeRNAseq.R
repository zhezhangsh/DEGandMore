# Normalize a read count matrix from RNA-seq
NormalizeRNAseq <- function(cnt, len=NA, round=FALSE, 
                            methods=c("NormTotalCount", "NormMedian", "NormQQ", "NormUpperQuantile", "NormRLE", 
                                      "NormTMM", "NormDESeq", "NormLinear", "NormLoess", "NormCyclicLoess", "NormFPKM", "NormTPM")) {
  
  normLog <- function(mthd, mtrx) {
    sz <- sum(rowMeans(mtrx+1)/10^6); 
    d  <- log2((mtrx+0.5)/sz); 
    
    d  <- do.call(mthd, list(mtrx=d)); 
    
    d  <- exp(d*log(2)); 
    d  <- d*sz-0.5;
    d[d<0] <- 0;
    d
  }; 
  
  log.transform <- function(ct) {
    c <- ct; 
    c[c>1] <- 1; 
    c <- ct[rowSums(c)==ncol(c), ]; 
    q <- apply(cbind(rowMeans(c), c), 2, function(x) quantile(log2(x), probs=seq(0, 1, 0.01)));
    
    nsd <- apply(q[, -1], 2, function(y) {
      x <- q[, 1]; 
      sapply(1:length(x), function(i) {
        l <- lm(y[-i] ~ x[-i]); 
        p <- predict(l, newdata=data.frame(y=y));
        f <- coefficients(l); 
        p <- f[[1]]+f[[2]]*x; 
        (p[i]-y[i])/sd(abs(y[-i]-p[-i])); 
      }); 
    }); 
  }
  
  norm <- list(Original=cnt); 
  
  if ("NormTotalCount" %in% methods) norm$Total_Count             <- NormTotalCount(cnt); 
  if ("NormMedian" %in% methods) norm$Median                      <- NormMedian(cnt); 
  if ("NormQQ" %in% methods) norm$Quantile_Quantile               <- NormQQ(cnt, ref='median'); 
  if ("NormUpperQuantile" %in% methods) norm$Upper_Quantile       <- NormUpperQuantile(cnt); 
  if ("NormTMM" %in% methods) norm$Trimmed_Mean                   <- NormTMM(cnt); 
  if ("NormRLE" %in% methods) norm$Relative_Log                   <- NormRLE(cnt); 
  if ("NormDESeq" %in% methods) norm$DESeq                        <- NormDESeq(cnt); 
  if ("NormFPKM" %in% methods & length(len)==nrow(cnt)) norm$FPKM <- NormFPKM(cnt, len); 
  if ("NormTPM" %in% methods & length(len)==nrow(cnt)) norm$TPM   <- NormTPM(cnt, len); 
  if ("NormLoess" %in% methods) norm$Linear                       <- normLog('NormLinear', cnt); 
  if ("NormLoess" %in% methods) norm$Loess                        <- normLog('NormLoess', cnt); 
  if ("NormCyclicLoess" %in% methods) norm$Cyclic_Loess           <- normLog('NormCyclicLoess', cnt); 

  if (round) lapply(norm, round); 
  
  norm; 
}