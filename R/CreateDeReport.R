CreateDeReport<-function(inputs, is.url=TRUE, dir.temp='.',
                         rmd="https://raw.githubusercontent.com/zhezhangsh/DEGandMore/master/examples/DeReport/DeReport.Rmd") {
  # inputs    A list of input variables, or path to the file
  # is.url    If the .Rmd template is an URL; a local file otherwise
  # rmd       URL or path to the .Rmd template file
  
  if (is.character(inputs)) inputs<-readRDS(inputs);
  
  lb<-strsplit("org.Hs.eg.db; RCurl; formatR; GenomeInfoDb; XVector; methods; bitops; utils; tools; grDevices; zlibbioc; digest; biclust; jsonlite; evaluate; DEGandMore; RSQLite; gage; lattice; png; pathview; graph; KEGGgraph; rstudioapi; DBI; Rgraphviz; curl; yaml; parallel; rJava; knitr; httr; stringr; htmlwidgets; xlsxjars; Biostrings; S4Vectors; graphics; gtools; datasets; stats; caTools; IRanges; DT; stats4; grid; base; Biobase; R6; snow; AnnotationDbi; XML; rmarkdown; gdata; isa2; magrittr; htmltools; gplots; awsomics; modeltools; MASS; BiocGenerics; RankProd; KEGGREST; flexclust; colorspace; xlsx; KernSmooth; stringi", '; ')[[1]];
  pkgs<-rownames(installed.packages());
  lb<-lb[!(lb %in% pkgs)];
  if (length(lb) > 0) stop("Required packages not installed: ", paste(lb, collapse='; '));
  
  if (!is.url) {
    if (!RCurl::url.exists(rmd)) stop("Rmd file not exists: ", rmd);
  } else {
    ln<-RCurl::getURL(rmd)[[1]];
    writeLines(ln, paste(dir.temp, 'DeReport.Rmd', sep='/'));
  }
  
  g1.ind<-inputs$indexes[[1]];
  g2.ind<-inputs$indexes[[2]];
  g1.name<-inputs$names[1];
  g2.name<-inputs$names[2];
  g1.name<-gsub('-', '_', g1.name);
  g2.name<-gsub('-', '_', g2.name);
  path0<-paste(g1.name, g2.name, sep='-vs-');
  path<-paste(inputs$output, path0, sep='/');
  
  rmarkdown::render(paste(dir.temp, 'DeReport.Rmd', sep='/'), output_file=paste(path, 'index.html', sep='/'));
  
  paste(path, 'index.html', sep='/');
}
