# Use the identify_outlier.Rmd template to create a report of gene clustering analysis
CreateOiReport<-function(yml) {
  # yml     The .ymal file or a yaml list defines the inputs and parameters of the analysis
  
  if (class(yml) == 'character') {
    if (!file.exists(yml)) stop('Input file ', yml, ' not found\n'); 
    yml <- yaml::yaml.load_file(yml);  
  }
  
  library(awsomics);
  library(gplots);
  library(knitr);
  library(rmarkdown); 
  
  if (!file.exists(yml$output)) dir.create(yml$output, recursive = TRUE)
  
  fn.temp<-paste(yml$output, 'identify_outlier.Rmd', sep='/'); 
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
  
  writeLines(yaml::as.yaml(yml), paste(yml$output, 'identify_outlier.yml', sep='/'));   
  
  zip(paste(yml$output, '.zip', sep=''), yml$output, "-r9X", zip='zip'); 
  
  list(index=fn.html, zip=paste(yml$output, '.zip', sep=''), error=errors);
}