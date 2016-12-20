# Combine the mean differences of a set of two-group comparisons
CombineMeans <- function(comps, sm = c('MD', 'SMD', 'ROM')) {
  require(meta); 
  
  sm <- toupper(sm[1]); 
  if (sm!='SMD' & sm!='ROM') sm <- 'MD';
  
  n1 <- sapply(comps, function(c) length(c[[2]][[1]])); 
  n2 <- sapply(comps, function(c) length(c[[2]][[2]])); 
  m1 <- lapply(comps, function(c) rowMeans(c[[1]][, c[[2]][[1]], drop=FALSE])); 
  m2 <- lapply(comps, function(c) rowMeans(c[[1]][, c[[2]][[2]], drop=FALSE])); 
  s1 <- lapply(comps, function(c) apply(c[[1]][, c[[2]][[1]], drop=FALSE], 1, sd)); 
  s2 <- lapply(comps, function(c) apply(c[[1]][, c[[2]][[2]], drop=FALSE], 1, sd)); 
  
  ip <- list(m1, m2, s1, s2); 
  id <- Reduce('union', lapply(m1, names)); 
  mx <- matrix(NA, nr=length(id), nc=length(comps), dimnames=list(id, names(comps)));
  ip <- lapply(ip, function(ip) {
    for (i in 1:length(comps)) mx[names(ip[[i]]), i] <- ip[[i]]; 
    mx; 
  }); 
  
  stat <- sapply(1:length(id), function(i) {
    m1 <- ip[[1]][i, ]; 
    m2 <- ip[[2]][i, ]; 
    s1 <- ip[[3]][i, ];
    s2 <- ip[[4]][i, ]; 

    mt <- metacont(n2, m2, s2, n1, m1, s1, sm=sm); 
    c(mt$k, mt$TE.random, mt$pval.random, mt$TE.fixed, mt$pval.fixed, mt$Q, mt$df.Q, mt$TE); 
  });
  stat <- t(stat); 
  
  qp <- pchisq(stat[, 6], stat[, 7], lower.tail = TRUE, log.p = TRUE);
  qp[stat[,7]==0 | is.na(qp)] <- 0;
  qp <- exp(qp); 
  
  m0 <- rowMeans(stat[, 8:ncol(stat), drop=FALSE], na.rm=TRUE); 
  
  stat <- cbind(stat[, 1:5], m0, exp(qp), stat[, 8:ncol(stat)]); 
  cnm <- paste(sm, names(comps), sep='_'); 
  colnames(stat) <- c('N', 'M_Random', 'P_Random', 'M_Fixed', 'P_Fixed', 'M_Mean', 'P_Heter', cnm);
  rownames(stat) <- id; 
  
  stat; 
}