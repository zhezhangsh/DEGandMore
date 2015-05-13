# Helper functions for preparing GSEA runs

########################################################################################################################
# Make the .gct file from a data matrix
PrepareGSEA<-function(mtrx, out, desc=rep('na', nrow(mtrx))) {
  # mtrx    The data matrix whose row names will be used as the <NAME> column in output
  # out     Path and name of output .gct file, with or without the .gct extension
  # desc    The values of the <Description> column in .gct file, which is optional; values will be 'na' if not given;
  
  if (!grepl('.gct$', out, ignore.case=TRUE)) out<-paste(out, '.gct', sep='');
  
  tbl<-data.frame(NAME=rownames(mtrx), Description=desc, mtrx, row.names=1:nrow(mtrx), stringsAsFactors=FALSE);
  writeLines(c("#1.2", 
               paste(nrow(mtrx), ncol(mtrx), sep='\t'),
               paste(colnames(tbl), collapse='\t')), 
             out);
  write.table(tbl, out, append=TRUE, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE);
  
  out;
}

########################################################################################################################

