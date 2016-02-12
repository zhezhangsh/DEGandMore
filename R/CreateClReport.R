# Use the ClReport.Rmd template to create a report of gene clustering analysis
CreateClReport<-function(yml) {
  # yml     The yaml file or an yaml list defines the inputs and parameters of the analysis

  if (class(yml) == 'character') {
    if (!exists('fn.yaml')) stop('Input file', yml, 'not found\n'); 
    yml <- yaml::yaml.load_file(fn.yaml);  
  }

  library(awsomics);
  library(gplots);
  library(knitr);
  library(rmarkdown); 
  
  
  if (!file.exists(yml$output)) dir.create(yml$output, recursive = TRUE)
  
  fn.temp<-paste(yml$output, 'ClReport.Rmd', sep='/'); 
  if (yml$input$remote) {
    if (!RCurl::url.exists(yml$input$template)) stop("Template Rmd file ', yml$input$template, ' not exists\n");
    writeLines(RCurl::getURL(yml$input$template)[[1]], fn.temp);
  } else {
    file.copy(yml$input$template, fn.temp); 
  }
  
  if (!file.exists(yml$output)) dir.create(yml$output, recursive = TRUE);
  
  fn.html<-paste(yml$output, 'index.html', sep='/'); 
  
  errors<-try(rmarkdown::render(fn.temp, output_format="html_document", output_file="index.html", output_dir=yml$output, 
                    quiet=TRUE, envir=new.env()), silent=TRUE);
  
  fn<-strsplit(fn.yaml, '/')[[1]];
  fn<-fn[length(fn)]; 
  
  file.copy(fn, paste(yml$output, fn, sep='/')); 
  zip(paste(yml$output, '.zip', sep=''), yml$output, "-rJ9X", zip='zip'); 
  
  list(index=fn.html, zip=paste(yml$output, '.zip', sep=''), error=errors);
}