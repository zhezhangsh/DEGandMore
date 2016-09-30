# Install dependency required to run the DE methods
installed <- rownames(installed.packages()); 

# Install required cran packages
pkg1 <- c('devtools', 'XML', 'readxl', 'shinythemes');
pkg  <- setdiff(pkg1, installed);
if (length(pkg) > 0) {
  install.packages(pkg); 
}

# install Bioconductor packages
pkg2 <- c("ABSSeq", "ALDEx2", "BADER", "baySeq", "BGmix", "bridge", "DEGseq", "DESeq2", "EBSeq", "edgeR", "edgeRun", 
          "limma", "LMGene", "NOISeq", "plgem", "PoissonSeq", "RankProd", "RBM", "samr", "sSeq", "tweeDEseq"); 
pkg  <- setdiff(pkg2, installed);
if (length(pkg) > 0) {
  source("http://bioconductor.org/biocLite.R");
  biocLite(pkg); 
}
# Install older version, the newest version doesn't support paired test
install_url("https://www.bioconductor.org/packages/3.2/bioc/src/contrib/ABSSeq_1.6.1.tar.gz");

# install GitHub packages (force to install to keep)
pkg3 <- c('RoCA', 'rchive', 'awsomics', 'DEGandMore')
require(devtools); 
install_github('zhezhangsh/RoCAR');
install_github('zhezhangsh/rchive');
install_github('zhezhangsh/awsomicsR');
install_github('zhezhangsh/DEGandMore');

# Check if all packages can be loaded
pkgs   <- c(pkg1, pkg2, pkg3); 
loaded <- sapply(pkgs, function(pkg) require(pkg, character.only = TRUE));
if(length(pkgs[!loaded]) == 0) cat("Congratulations! All required packages have been installed and loaded.\n") else
  cat('Oops!, Some packages are not installed/loaded: ', paste(pkgs[!loaded], collapse='; '), '\n'); 

