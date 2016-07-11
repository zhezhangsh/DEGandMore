# Normalize a read count matrix from RNA-seq
NormalizeRNAseq <- function(cnt, len=NA,  
                            methods=c("NormTotalCount", "NormMedian", "NormQQ", "NormUpperQuantile", "NormRLE", 
                                      "NormTMM", "NormDESeq", "NormLoess", "NormCyclicLoess", "NormFPKM", "NormTPM")) {
  
  normLog <- function(mthd, mtrx) {
    sz <- sum(rowMeans(mtrx+1)/10^6); 
    d  <- log2((mtrx+0.5)/sz); 
    d  <- do.call(mthd, list(mtrx=d)); 
    d  <- exp(d*log(2)); 
    d  <- d*sz-0.5;
    d[d<0] <- 0;
    d
  }; 
  
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
  if ("NormLoess" %in% methods) norm$Loess                        <- normLog('NormLoess', cnt); 
  if ("NormCyclicLoess" %in% methods) norm$Cyclic_Loess           <- normLog('NormCyclicLoess', cnt); 

  norm; 
}