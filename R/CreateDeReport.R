# Use this template to create a report of differential gene expression 
CreateDeReport<-function(yml) {
  # yml     The yaml file or an yaml list defines the inputs and parameters of the analysis
  
  library(awsomics);
  library(gplots);
  library(knitr);
  library(rmarkdown); 
  
  if (class(yml) == 'character') {
    if (!file.exists(yml)) stop('Input file', yml, 'not found\n'); 
    yml <- yaml::yaml.load_file(yml);  
  }

  if (!file.exists(yml$output)) dir.create(yml$output, recursive = TRUE);
  
  if (yml$input$remote) {
    if (!RCurl::url.exists(yml$input$template)) stop("Template Rmd file ', yml$input$template, ' not exists\n");
    writeLines(RCurl::getURL(yml$input$template)[[1]], './DeReport.Rmd');
  } else {
    file.copy(yml$input$template, './DeReport.Rmd'); 
  }
  
  errors<-try(rmarkdown::render('DeReport.Rmd', output_format="html_document", output_file="index.html", output_dir='.', 
                                quiet=TRUE, envir=new.env()), silent=TRUE);
  
  fn.html<-paste(yml$output, 'index.html', sep='/'); 

  paste(path, 'index.html', sep='/');
}
