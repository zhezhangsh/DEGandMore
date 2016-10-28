DeRNAseq <- function(ct, grps, paired = FALSE, mthds = 0, min.count = 6, num.cluster=1, just.stat = TRUE,
                     norm.count = c('DESeq', 'TMM', 'RLE', 'QQ', 'UpperQuantile', 'Median', 'TotalCount'),
                     norm.logged = c('Loess', 'VST', 'Rlog', 'QQ', 'UpperQuantile', 'Median'), force.norm = FALSE) {
  
  require(DEGandMore);
  data(DeMethodMeta);
  
  ######################################################################################################
  # Internal function
  normalizeCnt <- function(cnt, mthds) NormWrapper(cnt, paste('Norm', mthds, sep=''));
  normalizeLog <- function(cnt, mthds) {
    if (tolower(mthds[1]) %in% c('rlog', 'vst')) NormWrapper(cnt, paste('Norm', mthds, sep='')) else {
      cnt[cnt==0] <- 1/3;
      logged <- log2(cnt);
      NormWrapper(logged, mthd=paste('Norm', mthds, sep=''));
    }
  }
  ######################################################################################################
  
  ################################################################
  # Method names
  if (identical(mthds, 0)) mthds <- DeRNAseqMethods(0) else 
    if (identical(mthds, 1)) mthds <- DeRNAseqMethods(1) else 
      if (identical(mthds, 2)) mthds <- DeRNAseqMethods(2) else 
        if (identical(mthds, 3)) mthds <- DeRNAseqMethods(3) else 
          mthds <- mthds[mthds %in% DeRNAseqMethods()];
  if (length(mthds) == 0) mthds <- DeRNAseqMethods();
  ################################################################
  
  ################################################################################################################
  # Make sure the 2 groups to be compared are well-defined
  grps <- lapply(grps[1:2], function(s) s[s %in% colnames(ct) | (s>0 & s<=ncol(ct))]);
  if (length(grps[[1]])<2 | length(grps[[2]])<2) stop('No enough samples in one or both groups; requires 2.\n');
  if (is.null(names(grps))) names(grps) <- c('A', 'B');
  if (is.na(names(grps)[1])) names(grps)[1] <- 'A';
  if (is.na(names(grps)[2])) names(grps)[1] <- 'B';
  ct <- as.matrix(ct[, c(grps[[1]], grps[[2]])]);
  ################################################################################################################
  
  ###########################################################################################################
  # Only test genes with required number of total read count in all tested samples
  d0 <- ct[rowSums(ct) >= min.count, c(grps[[1]], grps[[2]]), drop=FALSE]; 
  if (nrow(d0) < 1) stop('No genes have total read count more than ', min.count, '; reduce the cutoff.\n'); 
  ###########################################################################################################
  
  #####################################################################################################################
  # Prepare normalized data if some methods require it
  norm <- DeMethodMeta[mthds, 'Normalization'];
  logd <- DeMethodMeta[mthds, 'Logged'];
  
  if (length(norm[norm=='No' & logd=='No'])>0 | force.norm) {
    norm1 <- paste('Norm', norm.count[1], sep=''); 
    if (!(norm1 %in% NormMethods())) stop('Normalization method not available: ', sub('Norm', '', norm.count), '\n');
    d1 <- NormWrapper(d0, norm1); 
  }
  if (length(norm[logd == 'Yes'])>0 | force.norm) {
    norm2 <- paste('Norm', norm.logged[1], sep=''); 
    if (!(norm2 %in% NormMethods())) stop('Normalization method not available: ', sub('Norm', '', norm.logged), '\n');
    if (norm.logged[1] %in% c('VST', 'Rlog')) {
      d2 <- NormWrapper(d0, norm2); 
    } else {
      d2 <- d0;
      d2[d2==0] <- 1/3;
      d2 <- log2(d2); 
      d2 <- NormWrapper(d2, norm2);       
    }
  }; 
  #####################################################################################################################
  
  ####################################################################
  # Prepare inputs
  d <- lapply(mthds, function(m) {
    lg <- DeMethodMeta[m, 'Logged']=='Yes';
    nm <- DeMethodMeta[m, 'Normalization']=='Yes'; 
    if (!lg & nm) d <- d0 else if (!lg) d <- d1 else d <- d2;
    list(method=m, data=d, group=grps, paired=paired, logged=lg); 
  });
  names(d) <- mthds;
  ####################################################################
  
  cat("Running DE analysis with methods: ", paste(mthds, collapse='; '), '\n'); 
  
  ###################################################################################################
  ###################################################################################################
  if (length(mthds)>1 & num.cluster>1 & require(snow)) {
    ####################################################################
    runDe <- function(d) {
      require(DEGandMore);
      DeWrapper(d$data, d$group, d$method, d$paired, d$logged)$results;
    }
    ####################################################################
    
    d  <- d[rev(order(DeMethodMeta[names(d), 'Speed']))];
    nm <- names(d); 
    
    cl <- snow::makeCluster(num.cluster, type='SOCK');
    stat <- snow::clusterApplyLB(cl, d[nm!='DePlgem'], runDe); 
    snow::stopCluster(cl); 
    names(stat) <- names(d)[nm!='DePlgem']; 
    
    if ('DePlgem' %in% nm) stat$DePlgem <- runDe(d[['DePlgem']]);
  } else {
    stat <- lapply(d, function(d) DeWrapper(d$data, d$group, d$method, d$paired, d$logged)$results);
  }
  stat <- stat[mthds]; 
  ###################################################################################################
  ###################################################################################################
  
  if (length(logd[logd=='Yes']) > 0) {
    un <- 2^d2;
    un <- un * (mean(d0, na.rm=TRUE)/mean(un, na.rm=TRUE));
    m1 <- rowMeans(un[, grps[[1]], drop=FALSE]);
    m2 <- rowMeans(un[, grps[[2]], drop=FALSE]);
    tb <- cbind(m1, m2, m2-m1); 
    stat[logd=='Yes'] <- lapply(stat[logd=='Yes'], function(s) {
      s$stat[, 1:3] <- tb;
      s;
    });
  }
  
  if (just.stat) stat <- lapply(stat, function(s) s$stat[, 1:6]); 
  
  normalized <- list();
  if (length(norm[norm=='No' & logd=='No'])>0 | force.norm) normalized$count  <- d1;
  if (length(logd[logd=='Yes'])>0 | force.norm) normalized$logged <- d2;
  
  input <- list(original=ct, filtered=d0, normalized=normalized, methods=mthds, groups=grps, paired=paired, 
                minimal.count=min.count, number.cluster=num.cluster, normalization=c(norm.count[1], norm.logged[1]));
  list(input=input, output=stat); 
}
##########################################################################################################

##########################################################################################################
# Method groups for DE analysis of RNA-seq data
DeRNAseqMethods <- function(group=NA) {
  
  require(DEGandMore);
  data(DeMethodMeta);
  
  ########################################################
  # Method categories
  default <- rownames(DeMethodMeta)[DeMethodMeta$Default=='Yes'];
  speed <- list(
    fast    = rownames(DeMethodMeta)[DeMethodMeta$Speed=='Fast'],
    medium  = rownames(DeMethodMeta)[DeMethodMeta$Speed=='Medium'],
    slow    = rownames(DeMethodMeta)[DeMethodMeta$Speed=='Slow'],
    slower  = rownames(DeMethodMeta)[DeMethodMeta$Speed=='Slower']
  );
  ########################################################
  
  group <- group[1];
  if (is.na(group)) group <- -1; 
  
  if (group==0) default else
    if (group==1) as.vector(unlist(speed[1])) else
      if (group==2) as.vector(unlist(speed[1:2])) else
        if (group==3) as.vector(unlist(speed[1:3])) else 
          unique(as.vector(unlist(speed)));
}
