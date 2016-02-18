# Helper functions for preparing GSEA runs

########################################################################################################################
# Make the .gct and .cls files from a data matrix
PrepareGSEA<-function(mtrx, grps, out, desc=rep('na', nrow(mtrx)), reverse.group=TRUE) {
  # mtrx    The data matrix whose row names will be used as the <NAME> column in output
  # grps    Named list of sample groups; each element has the sample indexes or names corresponding to columns in mtrx
  # out     Path and name prefix of output .gct file, with or without the .gct extension
  # desc    The values of the <Description> column in .gct file, which is optional; values will be 'na' if not given;

  mtrx<-mtrx[, unlist(grps, use.names=FALSE)];

  if (reverse.group) grps<-grps[2:1]; 
    
  # Prepare and write out the .gct file
  out1<-paste(out, '.gct', sep='');  
  tbl<-data.frame(NAME=rownames(mtrx), Description=desc, mtrx, row.names=1:nrow(mtrx), stringsAsFactors=FALSE);
  writeLines(c("#1.2", 
               paste(nrow(mtrx), ncol(mtrx), sep='\t'),
               paste(colnames(tbl), collapse='\t')), 
             out1);
  write.table(tbl, out1, append=TRUE, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE);

  # Prepare and write out the .cls file
  out2<-paste(out, '.cls', sep='');
  l1<-paste(ncol(mtrx), length(grps), 1, sep=' ');
  l2<-paste('#', paste(names(grps), collapse=' '), sep=' ');
  l3<-paste(rep(1:length(grps), sapply(grps, length))-1, collapse=' ');
  writeLines(c(l1, l2, l3), out2);
  
  c(out1, out2);
}

########################################################################################################################
