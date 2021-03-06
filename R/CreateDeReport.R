# Use this template to create a report of differential gene expression 
CreateDeReport<-function(yml, overwrite=FALSE) {
  # yml     The yaml file or an yaml list defines the inputs and parameters of the analysis
  
  library(awsomics);
  library(gplots);
  library(knitr);
  library(rmarkdown); 
  library(yaml);
  
  if (class(yml) == 'character') {
    if (!file.exists(yml)) stop('Input file', yml, 'not found\n'); 
    yml <- yaml::yaml.load_file(yml);  
  }

  if (!file.exists(yml$output)) dir.create(yml$output, recursive = TRUE);
  
  path<-paste(yml$output, paste(names(yml$input$comparison), collapse='-vs-'), sep='/');
  if (overwrite & file.exists(path)) unlink(path, recursive = TRUE); 
  
  if (yml$input$remote) {
    if (!RCurl::url.exists(yml$input$template)) stop("Template Rmd file ', yml$input$template, ' not exists\n");
    writeLines(RCurl::getURL(yml$input$template)[[1]], './DeReport.Rmd');
  } else {
    file.copy(yml$input$template, './DeReport.Rmd'); 
  }

  errors<-list(); 
  
  fn.md<-paste(path, 'DeReport.md', sep='/'); 
  if (file.exists(fn.md)) file.remove(fn.md); 
  
  # Run template, save error message
  errors$knit<-try(knit('DeReport.Rmd', fn.md)); 
  
  # Convert markdown file to html file
  errors$render<-try(rmarkdown::render(fn.md, output_format="html_document", output_file="index.html", output_dir=path, 
                                quiet=TRUE, envir=new.env()), silent=TRUE);
  
  # save template and yaml file
  file.copy('./DeReport.Rmd', paste(path, 'DeReport.Rmd', sep='/'));
  writeLines(as.yaml(yml), paste(path, 'DeReport.yml', sep='/'));
  
  # zip whole folder
  fn.zip<-paste(path, '.zip', sep=''); 
  zip(fn.zip, path, flag='-r0X', zip='zip');
  file.rename(fn.zip, paste(path, rev(strsplit(fn.zip, '/')[[1]])[1], sep='/')); 
  
  list(output=path, message=errors); 
}
