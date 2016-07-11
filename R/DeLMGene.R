# LMGene
# http://www.bioconductor.org/packages/release/bioc/html/LMGene.html
# LMGene {LMGene}
DeLMGene <- function(mtrx, grps, paired=FALSE, min.count=min(8, ncol(mtrx))) {
  
  require(DEGandMore);
  require(LMGene); 
  require(Biobase);
  require(tools); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];

  n        <- sapply(grps, length);

  if (paired & n[1]==n[2]) {
    prs <- rep(paste('pair', 1:n[1], sep='_'), 2);
    cnd <- rep(names(grps), n); 
    df  <- data.frame(condition=cnd, pair=prs); 
    mt  <- data.frame(labelDescription=c('condition', 'pair'));
  } else {
    cnd <- rep(names(grps), n); 
    df  <- data.frame(condition=cnd); 
    mt  <- data.frame(labelDescription=c('condition'));
  }
  lvls <- lapply(names(df), function(nm) df[[nm]]); 
  rownames(df) <- colnames(mtrx); 
  names(lvls) <- names(df); 
  
  anno <- AnnotatedDataFrame(data=df, varMetadata = mt); 
  eset <- ExpressionSet(mtrx, phenoData = anno); 
  esti <- tranest(eset, mult = TRUE, method=2);
  tran <- transeS(neweS(mtrx, lvls), lambda = esti[[1]], alpha = esti[[2]]); 

  res  <- genediff(tran, verbose = FALSE); 
  
  ex <- exprs(tran); 
  pv <- res$Gene.Specific[, 1];
  m1 <- rowMeans(ex[, grps[[1]], drop=FALSE]); 
  m2 <- rowMeans(ex[, grps[[2]], drop=FALSE]); 
  l2 <- m2-m1;
  
  pv[rowSums(mtrx) <= min.count] <- 1;
  qv <- p.adjust(pv, method='BH');
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps)[2:1], collapse='-'), 
                   'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, LMGene=res);
}