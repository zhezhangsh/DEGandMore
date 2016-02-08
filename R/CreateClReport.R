# Use the ClReport.Rmd template to create a report of gene clustering analysis
CreateClReport<-function(fn.yaml) {
  # fn.yaml     The .ymal file defines the inputs and parameters of the analysis

  if (!exists('fn.yaml')) stop('Input file not found\n'); 
  
  library(awsomics);
  library(gplots);
  library(knitr);
  library(rmarkdown); 
  
  yml <- yaml::yaml.load_file(fn.yaml);  
  
  fn.temp<-paste(yml$output, 'ClReport.Rmd', sep='/'); 
  if (yml$input$remote) {
    if (!RCurl::url.exists(yml$input$template)) stop("Template Rmd file ', yml$input$template, ' not exists\n");
    writeLines(RCurl::getURL(yml$input$template)[[1]], fn.temp);
  } else {
    file.copy(yml$input$template, fn.temp); 
  }
  
  fn.html<-paste(yml$output, 'index.html', sep='/'); 
  
  rmarkdown::render(fn.temp, output_file=fn.html);
  
  fn.html;
}