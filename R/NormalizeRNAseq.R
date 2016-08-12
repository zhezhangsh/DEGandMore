# Normalize a read count matrix from RNA-seq
NormalizeRNAseq <- function(cnt, len=NA, unlog=FALSE, round=FALSE, 
                            methods=c("NormCyclicLoess", "NormDESeq", "NormFPKM", "NormLinear", "NormLoess", "NormMedian", "NormQQ", 
                                      "NormRLE", "NormRlog", "NormTMM", "NormTotalCount", "NormTPM", "NormUpperQuantile", "NormVST")) {
  require(DEGandMore); 
  
  norm <- list(Original=cnt); 

  if ("NormCyclicLoess" %in% methods)                   norm$CL     <- NormCyclicLoess(log2(cnt+1)); 
  if ("NormDESeq" %in% methods)                         norm$DESeq  <- NormDESeq(cnt); 
  if ("NormFPKM" %in% methods & length(len)==nrow(cnt)) norm$FPKM   <- NormFPKM(cnt, len); 
  if ("NormLinear" %in% methods)                        norm$Linear <- NormLinear(log2(cnt+1)); 
  if ("NormLoess" %in% methods)                         norm$Loess  <- NormLoess(log2(cnt+1)); 
  if ("NormMedian" %in% methods)                        norm$Median <- NormMedian(cnt); 
  if ("NormQQ" %in% methods)                            norm$QQ     <- NormQQ(cnt); 
  if ("NormRLE" %in% methods)                           norm$RLE    <- NormRLE(cnt); 
  if ("NormRlog" %in% methods)                          norm$Rlog   <- NormRlog(cnt); 
  if ("NormTMM" %in% methods)                           norm$TMM    <- NormTMM(cnt); 
  if ("NormTotalCount" %in% methods)                    norm$TC     <- NormTotalCount(cnt); 
  if ("NormTPM" %in% methods & length(len)==nrow(cnt))  norm$TPM    <- NormTPM(cnt, len); 
  if ("NormUpperQuantile" %in% methods)                 norm$UQ     <- NormUpperQuantile(cnt); 
  if ("NormVST" %in% methods)                           norm$VST    <- NormVST(cnt); 
  
  if (unlog) {
    if (!is.na(norm$CL))     norm$CL     <- 2^norm$CL - 1;
    if (!is.na(norm$Linear)) norm$Linear <- 2^norm$Linear - 1;
    if (!is.na(norm$Loess))  norm$Loess  <- 2^norm$Loess - 1;
    if (!is.na(norm$Rlog))   norm$Rlog   <- 2^norm$Rlog;
    if (!is.na(norm$VST))    norm$VST    <- 2^norm$VST1;
    
    if (round) norm <- lapply(norm, round); 
  }
  
  norm; 
}



