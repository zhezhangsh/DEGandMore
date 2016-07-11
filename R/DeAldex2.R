# ALDEx2
# https://bioconductor.org/packages/release/bioc/html/ALDEx2.html
# aldex {ALDEx2}
DeAldex2 <- function(mtrx, grps, paired=FALSE, logged=TRUE, test=c('Welch', 'Wilcoxon'), mc.samples = min(128, 8*ncol(mtrx)), ...) {

  require(DEGandMore);
  require(ALDEx2);
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  norm <- aldex.clr(data.frame(mtrx), mc.samples=mc.samples, verbose=FALSE); 
  cond <- rep(names(grps), sapply(grps, length)); 
  
  res <- aldex.ttest(norm, cond, paired.test = paired); 
  res <- res[rownames(mtrx), ]; 
  
  eff <- aldex.effect(norm, cond, verbose = FALSE); 
  eff <- eff[rownames(mtrx), ]; 
  
  rownames(eff) <- rownames(res) <- rownames(mtrx); 
  
  if (tolower(test)[1] == 'wilcoxon') {
    pv <- res[, 3];
    p0 <- 'we.ep';
    p1 <- 'P_Welch';
  } else {
    pv <- res[, 1];
    p0 <- 'wi.ep';
    p1 <- 'P_Wilcoxon';
  }
  
  qv <- p.adjust(pv, method='BH');
  m1 <- eff[, 2];
  m2 <- eff[, 3]; 
  if (logged) l2 <- m2-m1 else {
    mn <- min(mtrx[mtrx>0]); 
    l2 <- log2(pmax(mn, m2)) - log2(pmax(mn, m1)); 
  };
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps)[2:1], collapse='-'), 
                   'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  s <- cbind(s, res[rownames(mtrx), p0], eff[, 'diff.btw'], eff[, 'effect']); 
  colnames(s)[7:9] <- c(p1, 'Median_Difference', 'Effect_Size'); 
  
  list(stat=s, group=grps, effect=eff);
}; 

