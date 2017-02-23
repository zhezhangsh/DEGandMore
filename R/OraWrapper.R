# Run over-representation test
OraWrapper<- function(subset, geneset, background=NA, method=c('fisher'=1, 'chisquare'=2, 'proportion'=3)) {
  
  gs1 <- unlist(geneset, use.names = FALSE);
  gs2 <- rep(names(geneset), sapply(geneset, length));
  
  if (identical(NA, background)) background <- unique(gs1);
  
  bgd <- background[!is.na(background) & background!=''];
  ind <- gs1 %in% bgd;
  gs1 <- gs1[ind];
  gs2 <- gs2[ind];
  gns <- subset[subset %in% bgd];
  
  if (length(gs1) < 1) c('No genes in gene sets found in background'=NULL) else
    if (length(gns) < 1) c('No genes in gene list found in background'=NULL) else {
      n11 <- n01 <- n10 <- n00 <- rep(0, length(geneset));
      names(n11) <- names(n01) <- names(n10) <- names(n00) <- names(geneset);
      in1 <- gs1 %in% gns;
      in0 <- !(gs1 %in% gns);
      mp1 <- split(gs1[in1], gs2[in1]);
      mp0 <- split(gs1[in0], gs2[in0]);
      n11[names(mp1)] <- sapply(mp1, length);
      n01[names(mp0)] <- sapply(mp0, length);
      n10[1:length(n10)] <- length(gns) - n11;
      n00[1:length(n00)] <- length(bgd) - length(gns) - n01;
      nnn <- cbind(n11, n01, n10, n00);
      
      osr <- n11*n00/(n10*n01);
      
      if (method == 1) {
        pvl <- phyper(n11-1, n11+n01, n10+n00, n11+n10, lower.tail=FALSE, log.p=TRUE);
        pvl <- exp(pvl); 
      } else if (method == 2) {
        options(warn = -1); 
        pvl <- apply(nnn, 1, function(ns) chisq.test(matrix(ns, nr=2))$p.value[[1]]);
        options(warn = -1); 
      } else if (method==3) {
        options(warn = -1); 
        pvl <- apply(nnn, 1, function(ns) prop.test(c(ns[1], ns[2]), c(ns[1]+ns[3], ns[2]+ns[4]), alternative='greater')$p.value[[1]])
        options(warn = 1); 
      } else pvl <- rep(1, length(pvl));
      
      tbl <- cbind(N_Within=n11, N_Total=n11+n01, "Within(%)"=round(100*n11/(n10+n10), 2), 
                   "Overall(%)"=round(100*(n11+n01)/rowSums(nnn), 2), "OddRatio"=round(osr, 2), 
                   PValue=pvl, 'FDR'=p.adjust(pvl, method='BH'));
      rownames(tbl) <- rownames(nnn); 
      
      tbl;
  }
  
  

}