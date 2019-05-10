# Create an index table for DE results generated from a number of groups
IndexDeRes <- function(groups, comps) {
  # groups:  names of all groups
  # comps:   a 3+ column data.frame. 1st column: name of group 0; 2nd column: name of group 1; and 3rd+ column: path to results
  
  cmp <- comps[comps[, 1] %in% groups & comps[, 2] %in% groups, , drop=FALSE]; 
  
  tbl <- matrix('', nrow = length(groups), ncol = length(groups), dimnames = list(groups, groups));
  
  if (nrow(cmp)>0) for (i in 1:nrow(cmp)) {
    cnm <- colnames(comps)[3:ncol(comps)]; 
    tbl[cmp[i, 1], cmp[i, 2]] <- paste(paste0('[', cnm, '](', cmp[i, 3:ncol(cmp)], '/index.html)'), collapse='|');
  }
  
  tbl; 
}