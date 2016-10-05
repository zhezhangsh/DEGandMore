# NOISeq
# https://www.bioconductor.org/packages/release/bioc/html/NOISeq.html
DeNoiSeq <- function(mtrx, grps, paired=FALSE, replicates=c('biological', 'technical', 'no'), 
                     norm = c("tmm", "n", "uqua", "rpkm"), ...) {
  # paired  Has no effect, paired test no supported
  
  require(DEGandMore);
  require(NOISeq);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];

  if (paired) warning("Paired test not supported by NOISeq; performing unpaired test instead.\n");
  
  rp <- tolower(replicates)[1]; 
  if (rp!='no' & rp!='technical') rp <- 'biological';
  
  nm <- tolower(norm)[1]; 
  if (nm!='n' & nm!='uqua' & nm!='rpkm') nm <- 'tmm';
  
  n <- sapply(grps, length); 
  f <- data.frame(comp=rep(names(grps), n)); 
  
  d <- readData(mtrx, factors=f, biotype = rp); 
  
  if (rp == 'biological') {
    noi <- noiseqbio(d, norm=nm, factor='comp', conditions = names(grps));
    res <- noi@results[[1]]; 
    
    m1 <- res[, 1];
    m2 <- res[, 2];
    l2 <- -res[, 5];
    pv <- 1-res[, 4];
    qv <- p.adjust(pv, method='BH');
  } else {
    noi <- noiseq(d, norm = nm, replicates = rp, factor="comp", conditions = names(grps));
    res <- noi@results[[1]]; 
    
    m1 <- res[, 1];
    m2 <- res[, 2];
    l2 <- log2(pmax(0.5, m2)) - log2(pmax(0.5, m1)); 
    pv <- 1-res[, 5];
    qv <- p.adjust(pv, method='BH');
  }; 
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(res); 
  
  list(stat=s[rownames(mtrx), ], group=grps, noi=noi); 
}