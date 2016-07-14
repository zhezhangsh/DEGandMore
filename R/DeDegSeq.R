# DEGSeq
# https://www.bioconductor.org/packages/release/bioc/html/DEGseq.html
# DEGexp {DEGseq}
DeDegSeq <- function(mtrx, grps, paired=FALSE, method=c("LRT", "CTR", "FET", "MARS", "MATR", "FC")) {
  
  require(DEGandMore);
  require(DEGseq);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by DEGseq; performing unpaired test instead.\n");
  
  mthd <- method[1]; 
  norm <- normalMethod[1];                     
  if (!(mthd %in% c("CTR", "FET", "MARS", "MATR", "FC"))) mthd <- 'LRT';      
  #if (!(norm %in% c("loess", "loess"))) norm <- 'none';

  mtrx1 <- as.matrix(cbind(rownames(mtrx), mtrx[, grps[[1]], drop=FALSE])); 
  mtrx2 <- as.matrix(cbind(rownames(mtrx), mtrx[, grps[[2]], drop=FALSE])); 
  
  DEGexp(geneExpMatrix1 = mtrx1, geneCol1=1, expCol1=2:ncol(mtrx1), groupLabel1 = names(grps)[1], 
         geneExpMatrix2 = mtrx2, geneCol2=1, expCol2=2:ncol(mtrx2), groupLabel2 = names(grps)[2],
         method = mthd, normalMethod = 'none', outputDir=tempdir()); 
  
  res <- read.table(paste(tempdir(), 'output_score.txt', sep='/'), row=1, header=TRUE, sep='\t'); 
  
  m1 <- res[, 1]/length(grps[[1]]);
  m2 <- res[, 2]/length(grps[[2]]);
  l2 <- res[, 4]; 
  pv <- res[, 5];
  qv <- p.adjust(pv, method='BH');
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps)[2:1], collapse='-'), 
                   'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(res); 
  
  lst <- list(geneExpMatrix1 = mtrx1, geneCol1=1, expCol1=2:ncol(mtrx1), groupLabel1 = names(grps)[1], 
              geneExpMatrix2 = mtrx2, geneCol2=1, expCol2=2:ncol(mtrx2), groupLabel2 = names(grps)[2],
              method = mthd, normalMethod = norm, outputDir=tempdir());
  
  list(stat=s[rownames(mtrx), ], group=grps, degseq=lst); 
}