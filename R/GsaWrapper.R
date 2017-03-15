GsaWrapper <- function(geneset, stat) {
  stat[1:length(stat)] <- as.numeric(stat);
  stat <- stat[!is.na(stat)]; 
}