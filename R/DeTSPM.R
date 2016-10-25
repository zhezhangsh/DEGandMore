DeTSPM <- function(mtrx, grps, paired=FALSE) {
  require(DEGandMore);
  require(GPseq); 
  
  prepared <- PrepareDe(mtrx, grps, paired);
  mtrx     <- prepared[[1]];
  grps     <- prepared[[2]];
  paired   <- prepared[[3]];
  
  if (paired) warning("Paired test not supported by TSPM; performing unpaired test instead.\n");
  
  #######################################################################
  # Downloaded from http://www.stat.purdue.edu/~doerge/software/TSPM.R  
  #######################################################################

  ##-------------------------------------------------------------------
  ## 	Name: TSPM.R
  ##  	R code for the paper by Paul L. Auer and R.W. Doerge:
  ## 	"A Two-Stage Poisson Model for Testing RNA-Seq Data"
  ## 	Date: February 2011
  ## 	Contact: Paul Auer 		plivermo@fhcrc.org
  ##		 R.W. Doerge 		doerge@purdue.edu
  
  
  ## 	Example:
  ## 	counts <- matrix(0, nrow=1000, ncol=10)
  ## 	for(i in 1:1000){
  ##		lambda <- rpois(n=1, lambda=10)
  ##		counts[i,] <- rpois(n=10, lambda=lambda)
  ## 	}
  ## 	x1 <- gl(n=2, k=5, labels=c("T", "C"))
  ## 	x0 <- rep(1, times=10)
  ## 	lib.size <- apply(counts,2,sum)
  ##	result <- TSPM(counts, x1, x0, lib.size)
  ##---------------------------------------------------------------------
  
  
  #######################################################################
  ###### The TSPM function ##############################################
  #######################################################################
  TSPM <- function(counts, x1, x0, lib.size, alpha.wh=0.05){
    
    ## Input:
    #counts: 	a matrix of RNA-Seq gene counts (genes are rows, samples are columns)
    #x1: 		a vector of treatment group factors (under the alternative hypothesis)
    #x0: 		a vector of treatment group factors (under the null hypothesis)
    #lib.size: 	a vector of RNA-Seq library sizes. This could simply be obtained
    #          	by specifying lib.size <- apply(counts,2,sum). It may also be any other
    #          	appropriate scaling factor.
    #alpha.wh:	the significance threshold to use for deciding whether a gene is overdispersed.
    #               Defaults to 0.05.
    
    
    ## Output:
    #log.fold.change:		a vector containing the estimated log fold changes for each gene
    #pvalues: 			a vector containing the raw p-values testing differential expression for each gene.
    #index.over.disp: 	a vector of integer values containing the indices of the over-dispersed genes.
    #index.not.over.disp:	a vector of integer values containing the indices of the non-over-dispersed genes.
    #padj:			a vector containing the p-values after adjusting for multiple testing using the 
    #				method of Benjamini-Hochberg
    
    
    ######## The main loop that fits the GLMs to each gene ########################
    
    ### Initializing model parameters ####
    n <- dim(counts)[1]
    per.gene.disp <- NULL
    LRT <- NULL
    score.test <- NULL
    LFC <- NULL
    
    ###### Fitting the GLMs for each gene #################
    for(i in 1:n){
      ### Fit full and reduced models ###
      model.1 <- glm(as.numeric(counts[i,]) ~ x1, offset=log(lib.size), family=poisson)
      model.0 <- glm(as.numeric(counts[i,]) ~ x0, offset=log(lib.size), family=poisson)
      
      ### Obtain diagonals of Hat matrix from the full model fit ###
      hats <- hatvalues(model.1)
      
      ### Obtain Pearson overdispersion estimate ####
      per.gene.disp[i] <- sum(residuals(model.1, type="pearson")^2)/model.1$df.residual
      
      ### Obtain Likelihood ratio statistic ####
      LRT[i] <- deviance(model.0)-deviance(model.1)
      
      ### Obtain score test statistic ####
      score.test[i] <- 1/(2*length(counts[i,])) * sum(residuals(model.1, type="pearson")^2 - ((counts[i,] - hats*model.1$fitted.values)/model.1$fitted.values))^2
      
      ### Obtain the estimated log fold change ###
      LFC[i] <- -model.1$coef[2]
    }
    
    ## Initialize parameters for Working-Hotelling bands around the score TSs ###
    qchi <- qchisq(df=1, (1:n-0.5)/n)
    MSE <- 2
    UL <- NULL
    
    #### Obtain the upper boundary of the WH bands #######################################
    xbar <- mean(qchi)
    bottom <- sum((qchi-xbar)^2)
    top <- (qchi-xbar)^2
    s <- sqrt(MSE*(1/n) + (top/bottom))
    W <- sqrt(2*qf(df1=1, df2=n-1, p=1-(alpha.wh/n)))
    UL <- pmax(qchi + W*s,1)
    
    ###### Obtain the indices of the over-dispersed and not-over-dispersed genes, respectively ##########
    
    cutoff <- min(which(sort(score.test)-UL > 0))
    temp <- cutoff-1 + seq(cutoff:length(score.test))
    over.disp <- which(score.test %in% sort(score.test)[temp])
    not.over.disp <- setdiff(1:length(score.test), over.disp)
    
    ###### Compute p-values ####################################
    p.f <- pf(LRT[over.disp]/per.gene.disp[over.disp], df1=1, df2=model.1$df.residual, lower.tail=FALSE)
    p.chi <- pchisq(LRT[not.over.disp], df=1, lower.tail=FALSE)
    p <- NULL
    p[over.disp] <- p.f
    p[not.over.disp] <- p.chi
    
    ##### Adjust the p-values using the B-H method ####################
    p.bh.f <- p.adjust(p.f, method="BH")
    p.bh.chi <- p.adjust(p.chi, method="BH")
    final.p.bh.tagwise <- NULL
    final.p.bh.tagwise[over.disp] <- p.bh.f
    final.p.bh.tagwise[not.over.disp] <- p.bh.chi
    
    ### Output ###
    list(log.fold.change=LFC, pvalues=p, index.over.disp=over.disp, index.not.over.disp=not.over.disp,
         padj=final.p.bh.tagwise)
  }
  
  n  <- sapply(grps, length); 
  x1 <- as.factor(rep(names(grps), n));
  x0 <- rep(1, sum(n)); 
  
  ##########################################
  res <- TSPM(mtrx, x1, x0, colSums(mtrx)); 
  ##########################################
  
  l2 <- -res$log.fold.change;
  l2 <- log2(exp(l2)); 
  pv <- res$pvalues;
  qv <- p.adjust(pv, method='BH');
  
  aj <- colSums(mtrx)/mean(colSums(mtrx));
  m1 <- rowMeans(mtrx[, grps[[1]]])/mean(aj[grps[[1]]]);
  m2 <- rowMeans(mtrx[, grps[[2]]])/mean(aj[grps[[2]]]);
  
  m1[is.na(m1)] <- 0;
  m2[is.na(m2)] <- 0;
  l2[is.na(l2)] <- 0;
  pv[is.na(pv)] <- 1;
  qv[is.na(qv)] <- 1;
  
  s <- cbind(m1, m2, m2-m1, l2, pv, qv);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), 'Mean_Change', 'LogFC', 'Pvalue', 'FDR');
  rownames(s) <- rownames(mtrx); 
  
  list(stat=s, group=grps, TSPM=res);
}