# Functions used for making PCA plots of compared sample groups

# Generic plotting of PCA results
DegPCA<-function(pca, groups, filename=NA, dimensions=1:2, legend=FALSE, legend.single=FALSE, 
                 col=c(), pch=19, cex=2, labels=NA, col.labels='white', propagation=FALSE, new.window=TRUE) {
  
  # pca			  Outputs of prcomp() or a data matrix, whose columns are sample IDs
  # groups		The grouping label of the samples
  # filename		Output to a pdf file if not NA
  # groupnames	The names of the 2 groups of samples
  # legend=FALSE	Whether to draw a legend; if TRUE, split the canvas and use the second panel and function PCALegend to draw legend
  # dimensions	The 2 dimensions of PCA outputs to be used
  # col,pch,cex	Plot parameters
  # labels		If have the same length as the number of samples, label the samples
  # col.labels	Color of the sample labels
  # propagation   Whether to run an affinity propagation clusterings
  # new.window    If TRUE, plot the PCA result in a new Quartz
  #
  # RETURN		A list of 3 elements
  #				[[1]] the current plotting device
  #				[[2]] a data frame, each row for a sample, with these fields in order: group label, col, pch, cex, and sample labels if individual samples were labelled
  #				[[3]] a data frame, each row for a group, with these fields in order: group name, col, pch, and cex
  
  CL<-c('black', 'red', 'green3', 'blue', 'cyan', 'magenta', 'orange', 'gray', 'antiquewhite', 'antiquewhite4', 'aquamarine', 'darkgreen', 'darkblue', 
        'blueviolet', 'brown', 'chartreuse1', 'yellow', 'darkseagreen1', 'deeppink1', 'gold', 'goldenrod4', 'darkmagenta', 'paleturquoise1', 'indianred1', 
        'hotpink', 'purple4', 'grey39', 'orchid4', 'turquoise4', 'yellow3', 'royalblue4', 'darkorange');
  
  if (class(pca) != 'prcomp') pca<-prcomp(t(pca));
  
  X<-pca$x[,dimensions[1]];
  Y<-pca$x[,dimensions[2]];
  per<-round(summary(pca)$importance[2, dimensions]*100, 2);
  
  groups<-as.vector(groups);
  groupnames<-unique(groups)
  
  
  # set col
  if (length(col)!=length(X)) {
    if (length(col)==1) col<-rep(col, length(X))
    else if (length(col)==length(groupnames)) {
      names(col)<-groupnames;
      col<-col[groups];
    }
    else {
      col<-CL[1:length(groupnames)];
      names(col)<-groupnames;
      col<-col[groups];
    }
  } # end of col
  
  # set pch
  if (length(pch)!=length(X)) { # if setting for each sample not available
    if (length(pch)==1) pch=rep(pch, length(X))
    else if (length(pch)==length(groupnames)) {
      names(pch)<-groupnames;
      pch<-pch[groups];
    }
    else pch<-rep(19, length(X));
  } # end of pch
  
  # set cex
  if (length(cex)!=length(X)) { # if setting for each sample not available
    if (length(cex)==1) cex=rep(cex, length(X))
    else if (length(cex)==length(groupnames)) {
      names(cex)<-groupnames;
      cex<-cex[groups];
    }
    else cex<-rep(2, length(X));
  } # end of cex
  
  names(col)<-names(cex)<-names(pch)<-names(X);
  
  # set output device
  if (legend) w=10 else w=7.5;
  h=7.5;
  if (!identical(NA, filename)) {
    if (gregexpr('pdf$', filename, ignore.case=TRUE)==-1) filename=paste(filename, '.pdf', sep='');
    pdf(filename, w=w, h=h) # write to pdf file
  }
  else if (new.window) quartz(w=w, h=h);
  
  if (legend) layout(matrix(1:2, nrow=1), width=c(3,1));
  par(mai=c(1,1,.25,.25));
  
  plot(X, Y, col=col, pch=pch, cex=cex, xlim=c(min(X)*1.1, max(X)*1.1), ylim=c(min(Y)*1.1, max(Y)*1.1), 
       xlab=paste('PC', dimensions[1], ', ', per[1], '%', sep=''), ylab=paste('PC', dimensions[2], ', ', 
                                                                              per[2], '%', sep=''), cex.lab=2.5);
  
  if (propagation) {
    AddAPCluster(cbind(X,Y), cex);
  }
  
  # draw labels
  if (length(labels)==length(X)) {
    text(X, Y, label=labels, col=col.labels, cex=cex/3);
  }
  
  x<-data.frame(id=groups, col=col, pch=pch, cex=cex);
  if (length(labels)==length(X)) x<-cbind(x, labels);
  
  #if (legend) cat('Please use information in outputs and function PCALegend to draw legend')
  if (legend) {
    if (legend.single) PCALegend(dev.cur(), settings=x)
    else PCALegend(dev.cur(), settings=unique(x[, 1:4]));
  }
  else if (!identical(NA, filename))  capture.output(z<-lapply(dev.list(), function(x) dev.off(which=x)));
  
  filename;
}
#################################################################################


#################################################################################
# A simplfied version of DegPCA, which draw PCA of 2 groups of samples
DegPCA2Groups<-function(pca, ind1, ind2, label.samples=FALSE, filename=NA, groupnames=c('A', 'B'), 
                         dimensions=1:2, col=c(4, 2), pch=c(19, 19), cex=c(2,2), propagation=FALSE, new.window=TRUE) {
  # pca  		    Outputs of prcomp() or a data matrix, whose columns are sample IDs
  # ind1, ind2	Indexes of samples of the 2 groups
  # label.samples	Whether to label indivdual samples with ordered integer; if FALSE, only put group names in legend
  # filename		Output to a pdf file if not NA
  # groupnames	The names of the 2 groups of samples
  # dimensions	The 2 dimensions of PCA outputs to be used
  # col,pch,cex	Plot parameters
  # new.window    If TRUE, plot the PCA result in a new Quartz

  ind<-c(ind1, ind2);

  if (class(pca) != 'prcomp') pca<-prcomp(t(pca[, ind]));
  
  X<-pca$x[ind, dimensions[1]];
  Y<-pca$x[ind, dimensions[2]];
  per<-round(summary(pca)$importance[2, dimensions]*100, 2);
  
  # set plot parameters
  col<-rep(col, c(length(ind1), length(ind2)));
  pch<-rep(pch, c(length(ind1), length(ind2)));
  cex<-rep(cex, c(length(ind1), length(ind2)));
  names(col)<-names(cex)<-names(pch)<-names(X);
  
  # set output device
  w=10
  h=7.5;
  if (!identical(NA, filename)) {
    if (gregexpr('pdf$', filename, ignore.case=TRUE)==-1) filename=paste(filename, '.pdf', sep='');
    pdf(filename, w=w, h=h) # write to pdf file
  }
  else if (new.window) quartz(w=w, h=h);
  layout(matrix(1:2, nrow=1), width=c(3,1));
  par(mai=c(1,1,.25,.25));
  
  
  plot(X, Y, col=col, pch=pch, cex=cex, xlim=c(min(X)*1.1, max(X)*1.1), ylim=c(min(Y)*1.1, max(Y)*1.1), 
       xlab=paste('PC', dimensions[1], ', ', per[1], '%', sep=''), ylab=paste('PC', dimensions[2], ', ', 
                                                                              per[2], '%', sep=''), cex.lab=2.5);
  
  if (propagation) {
    AddAPCluster(cbind(X,Y), cex);
  }
  
  if (label.samples) {
    l<-c(1:length(ind1), 1:length(ind2));
    text(X, Y, label=l, col='white', cex=cex/3);
    PCALegend(device=dev.cur(), data.frame(names(X), col, pch, cex, l));
  }
  else  {
    i<-c(1, length(ind1)+1);
    PCALegend(device=dev.cur(), data.frame(groupnames, col[i], pch[i], cex[i]));
  }
}

#################################################################################
#################################################################################
AddAPCluster<-function(d, cex=1) {
  #d   matrix of numbers, the first 2 columns will be used for plotting
  library(apcluster);
  
  if (ncol(d)<2) {
    cat('Only enough columns in input data for AP clustering\n');
  }
  else {
    s<-negDistMat(d, r = 2)
    ap<-apcluster(s); 
    e<-ap@exemplars};
  c<-ap@clusters;
  
  points(d[e,1], d[e,2], pch=22, cex=cex+1);
  
  for (i in 1:length(e)) {
    x0<-d[e[i],1];
    y0<-d[e[i],2];
    x1<-d[c[[i]],1];
    y1<-d[c[[i]],2];
    segments(x0, y0, x1, y1, lty=3, col='darkgrey');
  }
}



#################################################################################
#################################################################################
# plot the legend of a PCA
PCALegend<-function(device, settings, labels=as.vector(settings[,1])) {
  # device		The output device
  # settings		Data.frame, the contents of the legend, with the following fields in order: legend label, color, pch, size, and individual label of samples (optional)
  # labels		String vector, legend labels if different from the first field of settings
  
  dev.set(device);
  par(mai=c(1, 0, 0.25, 0));
  plot(0, type='n', xlim=c(0, 100), ylim=c(1, 100), axes=FALSE, bty='n', xaxs='i', yaxs='i', xlab='', ylab='');
  
  w<-strwidth(labels, cex=1.2);
  w0<-1.2*80/max(w, 80);
  h0<-1.2*(100/length(labels))/4.0
  cex<-min(1.2, w0, h0); 
  cex<-cex*settings[,4]/max(settings[,4]);
  
  points(rep(5, length(labels)), 100-(1:length(labels))*4.0*cex/1.25, col=as.vector(settings[,2]), pch=settings[,3], cex=1.5*cex);
  text(5+cex*5, 100-(1:length(labels))*4.0*cex/1.25, labels=as.vector(settings[,1]), adj=0, col=as.vector(settings[,2]), pch=settings[,3], cex=cex);
  if (ncol(settings)>4) text(rep(5, length(labels)), 100-(1:length(labels))*4.0*cex/1.25, labels=as.vector(settings[,5]), col='white', cex=cex/1.5);
  
  #legend(0, 100, legend=labels, bty='n', col=as.vector(settings[,2]), text.col=as.vector(settings[,2]), pch=as.vector(settings[,3]), pt.cex=as.vector(settings[,4])*min(1.5, cex), cex=cex);
  if(names(device)!='quartz') capture.output(x<-dev.off(which=device));
}