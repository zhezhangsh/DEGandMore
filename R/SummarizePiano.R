SummarizePiano <- function(s) {
  n <- sapply(s, length); 
  s <- s[n==n['gsc']];
  i <- names(s$gsc); 
  s <- s[names(s) != 'gsc']; 
  n <- sapply(s, function(s) length(s[!is.na(s)])); 
  n <- n[n==max(n)]; 
  o <- do.call('cbind', s[names(n)]); 
  colnames(o) <- names(n); 
  rownames(o) <- i; 
  o; 
}