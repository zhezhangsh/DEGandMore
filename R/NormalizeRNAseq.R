# Normalize a read count matrix from RNA-seq
NormalizeRNAseq <- function(cnt, len=NA, round=FALSE, 
                            methods=c("NormCyclicLoess", "NormDESeq", "NormFPKM", "NormLinear", "NormLoess", "NormMedian", "NormQQ", 
                                      "NormRLE", "NormRlog", "NormTMM", "NormTotalCount", "NormTPM", "NormUpperQuantile", "NormVST")) {
  require(DEGandMore); 
  
  # Perform normalizaiton on log-scale data matrix
  normLog <- function(mthd, mtrx) {
    d <- log2(mtrx+1); 
    d <- do.call(mthd, list(mtrx=d)); 
    d <- 2^d - 1; 
    d[d<0] <- 0;
    d
  }; 
  
  norm <- list(Original=cnt); 

  if ("NormCyclicLoess" %in% methods)                   norm$CL     <- normLog('NormCyclicLoess', cnt); 
  if ("NormDESeq" %in% methods)                         norm$DESeq  <- NormDESeq(cnt); 
  if ("NormFPKM" %in% methods & length(len)==nrow(cnt)) norm$FPKM   <- NormFPKM(cnt, len); 
  if ("NormLinear" %in% methods)                        norm$Linear <- normLog('NormLinear', cnt); 
  if ("NormLoess" %in% methods)                         norm$Loess  <- normLog('NormLoess', cnt); 
  if ("NormMedian" %in% methods)                        norm$Median <- NormMedian(cnt); 
  if ("NormQQ" %in% methods)                            norm$QQ     <- NormQQ(cnt); 
  if ("NormRLE" %in% methods)                           norm$RLE    <- NormRLE(cnt); 
  if ("NormRlog" %in% methods)                          norm$Rlog   <- 2^NormRlog(cnt); 
  if ("NormTMM" %in% methods)                           norm$TMM    <- NormTMM(cnt); 
  if ("NormTotalCount" %in% methods)                    norm$TC     <- NormTotalCount(cnt); 
  if ("NormTPM" %in% methods & length(len)==nrow(cnt))  norm$TPM    <- NormTPM(cnt, len); 
  if ("NormUpperQuantile" %in% methods)                 norm$UQ     <- NormUpperQuantile(cnt); 
  if ("NormVST" %in% methods)                           norm$VST    <- 2^NormVST(cnt); 
  
  if (round) norm <- lapply(norm, round); 
  
  norm; 
}