DeVoomLimma<-function(mtrx, grps, paired=FALSE, plot=FALSE, ...) {
  require(DEGandMore);
  require("edgeR");
  require("limma");

  e1   <- mtrx[, grps[[1]], drop=FALSE];
  e2   <- mtrx[, grps[[2]], drop=FALSE];
  mtrx <- cbind(e1, e2);
  m1   <- rowMeans(e1, na.rm=TRUE);
  m2   <- rowMeans(e2, na.rm=TRUE);
  
  res  <- VoomLimma(mtrx=mtrx, grps=grps, paired=paired, plot=plot, ...);
  
  lgfc <- res[[1]][, 1];
  p    <- res[[1]][, 2];
  q    <- res[[2]][, 'adj.P.Val'];
  lgfc[is.na(lgfc)] <- 0;
  p[is.na(p)]       <- 1;
  q[is.na(q)]       <- 1;
  
  s <- cbind(cbind(m1, m2, m2-m1)[rownames(res[[1]]), ], lgfc, p, q);
  colnames(s) <- c(paste('Mean', names(grps), sep='_'), paste(names(grps), collapse='-'), 'LogFC', 'Pvalue', 'FDR');
  s <- cbind(s, res[[2]][, c(2, 3, 6)]);
  
  list(stat=s[rownames(mtrx), ], group=grps, voom=res);
}

###########################################
#Code to run Voom and Limma on a htseq count matrix
#
#Function requires a matrix of counts and 2 or more groups
#returns back a list
#
#First object is gene, p-value, Log-FC
#Second object is the actual limma output
#
###########################################

#####################
#EXAMPLES
#
#2 Groups:
#load("../data/mouse.brain.htseq.rda");
#groupsName <- list(colnames(mouse.brain.htseq)[1:5], colnames(mouse.brain.htseq)[6:10]);
#groudpsNum <- list(c(1:5), c(6:10));
#limmaOutput2<- voomLimma(mouse.brain.htseq, groudpsNum);
#
#Multiple Groups:
#groupsName <- list(colnames(mouse.brain.htseq)[1:3], colnames(mouse.brain.htseq)[4:6], colnames(mouse.brain.htseq)[7:10]);
#limmaOutputMult<- voomLimma(mouse.brain.htseq, groupsName);
#
#####################

#Function to run voom and then limma, by default first group is control
#If more than one group is given the ANOVA p-value is outputter
#all other groups are always compared back to control or you get an explosion in number of comparisons
VoomLimma <- function(mtrx, grps, paired=FALSE, plot=FALSE, ...) {
  #Call libraries
  require("edgeR");
  require("limma");
  
  if (paired) if (length(grps[[1]]) != length(grps[[2]])) paired <- FALSE;
  
  dataMat <- mtrx;
  groups <- grps;
  #Handle if its a numeric vector
  if(class(groups[[1]])=="integer")
  {
    tmpGroups <- colnames(dataMat)[groups[[1]]];
    for(i in 2:length(groups))
    {
      tmpGroups <- list(tmpGroups, colnames(dataMat)[groups[[i]]]);
    }
    groups <- tmpGroups;
  }
  
  #Generate targets file
  targets <- data.frame(groups[[1]],paste("Group_", rep(1, length(groups[[1]])), sep=""))
  colnames(targets) <- c("Sample", "Condition");
  for(i in 2:length(groups))
  {
    tmpTargets <- data.frame(groups[[i]],paste("Group_", rep(i, length(groups[[i]])), sep=""))
    colnames(tmpTargets) <- c("Sample", "Condition");
    targets <- rbind(targets, tmpTargets);
  }
  rownames(targets) <- targets[,1];
  if (paired) targets$Pair <- paste('Pair', rep(1:(nrow(targets)/2), 2), sep='_')
  targets <- targets[, -1];
  
  #Run Voom
  y <- DGEList(counts=as.matrix(dataMat), genes=rownames(dataMat))
  y <- calcNormFactors(y)
  if (paired) design <- model.matrix(~Pair+Condition, data=targets) else 
    design <- model.matrix(~Condition, data=targets)
  v <- voom(y,design,plot=plot, ...);
  voomData <- v$E
  
  #Run linear fit and Create the contrast
#   if (paired) design <- model.matrix(~0+Condition+Pair, data=targets) else 
#     design <- model.matrix(~0+Condition, data=targets)
  design <- model.matrix(~0+targets[,"Condition"])
  colnames(design)[1:2] <- levels(targets[,1]);
  fit <- lmFit(voomData, design);
  myConts <- paste(levels(targets[,1])[-1], "Group_1", sep="-");
  contrast.matrix <- makeContrasts(contrasts=myConts, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  #Okay now combine all into one result data frame
  myComps <- colnames(fit2$contrasts);
  
  if(length(myComps)==1)
  {
    output <- topTable(fit2, number=nrow(mtrx))
    genOutput <- output[,c("logFC", "P.Value")];
    AllResult <- list(genOutput, output, voomData, design, myComps);
  }
  
  if(length(myComps)>1)
  {
    output <- topTable(fit2, number=nrow(mtrx))
    genOutput <- output[,c(c(1:length(myComps)), grep("P.Value", colnames(output)))];
    colnames(genOutput)[1:length(myComps)] <- myComps;
    tmpOut <- topTable(fit2, number=nrow(mtrx), 1)[,c("logFC", "P.Value", "adj.P.Val")];
    colnames(tmpOut) <- paste(myComps[1], colnames(tmpOut), sep="_");
    IndComp <- tmpOut;
    
    for(i in 2:length(myComps))
    {
      tmpOut <- topTable(fit2, number=nrow(mtrx), i)[,c("logFC", "P.Value", "adj.P.Val")];
      tmpOut <- tmpOut[rownames(IndComp),];
      colnames(tmpOut) <- paste(myComps[i], colnames(tmpOut), sep="_");
      IndComp <- cbind(IndComp, tmpOut);
    }
    
    AllResult <- list(genOutput, IndComp, voomData, design, myComps);
  }
  names(AllResult) <- c("effectPval", "fullResult", "voomMatrix", "design", "comparison")
  return(AllResult);
}