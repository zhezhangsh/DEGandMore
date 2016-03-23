# Use this template to create a report that compares 2 pairwise comparisons of differential gene expression 
CreateDdReport<-function(yml, overwrite=FALSE) {
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
  
  path<-yml$output;
  if (!file.exists(path)) dir.create(, recursive = TRUE);

  if (yml$template$remote) {
    if (!RCurl::url.exists(yml$template$location)) stop("Template Rmd file ', yml$input$template, ' not exists\n");
    writeLines(RCurl::getURL(yml$template$location)[[1]], './DdReport.Rmd');
  } else {
    file.copy(yml$template$location, './DdReport.Rmd'); 
  }

  errors<-list(); 
  
  fn.md<-paste(path, 'DdReport.md', sep='/'); 
  if (file.exists(fn.md)) file.remove(fn.md); 
  
  # Run template, save error message
  errors$knit<-try(knit('DdReport.Rmd', fn.md)); 
  
  file.rename('figure', paste(path, 'figure', sep='/')); 
  
  # Convert markdown file to html file
  errors$render<-try(rmarkdown::render(fn.md, output_format="html_document", output_file="index.html", output_dir=path, 
                                quiet=TRUE, envir=new.env()), silent=TRUE);
  
  # save template and yaml file
  file.copy('./DdReport.Rmd', paste(path, 'DdReport.Rmd', sep='/'));
  writeLines(as.yaml(yml), paste(path, 'DdReport.yml', sep='/'));
  
  # zip whole folder
  fn.zip<-paste(path, '.zip', sep=''); 
  zip(fn.zip, path, flag='-r0X', zip='zip');
  file.rename(fn.zip, paste(path, rev(strsplit(fn.zip, '/')[[1]])[1], sep='/')); 
  
  list(output=path, message=errors); 
}
