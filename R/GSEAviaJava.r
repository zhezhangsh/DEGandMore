# Run GSEA by creating a Java shell command line
GSEAviaJava<-function(yml, execute=TRUE) {
  # yml      yml file specifying all inputs or a corresponding yml list
  # execute   Whether to make the GSEA run or just save the command lines
  
  library(DEGandMore);
  library(yaml);
  
  if (class(yml) == 'character') {
    if (!file.exists(yml)) stop('Input file', yml, 'not found\n'); 
    yml<-yaml.load_file(yml);
  }
  prefix<-yml$name;
  if (!file.exists(prefix)) dir.create(prefix, recursive = TRUE);
  
  # build shell command line
  cmmd<-paste('java -Xmx5G -cp', yml$jar);
  
  # output folder
  out<-yml$out; # folder of run outputs
  if (is.null(out)) out<-prefix
  if (!file.exists(out)) dir.create(out, recursive=TRUE);
  
  # Copy input files to output folder
  fn0<-unlist(yml[c('input', 'class', 'chip', 'gmt', 'jar')]);
  fn1<-sapply(strsplit(fn0, '/'), function(fn) fn[length(fn)]);
  fn1<-paste(prefix, fn1, sep='/'); 
  names(fn1)<-names(fn0); 
  file.copy(fn0[fn0!=fn1], fn1[fn0!=fn1], overwrite = TRUE);
  
  # specify inputs
  if (yml$preranked) {
    cmmd<-paste(cmmd, 'xtools.gsea.GseaPreranked', '-rnk', fn1['input']);
  } else {
    cmmd<-paste(cmmd, 'xtools.gsea.Gsea', '-res', fn1['input'], '-cls', fn1['class']);
  }
  
  # add chip annotation
  if (!is.null(yml$chip)) if (file.exists(yml$chip)) cmmd<-paste(cmmd, '-chip', fn1['chip']);
  
  # extra options
  options<-yml$options;
  cmmd<-paste(cmmd, paste(sapply(names(options), function(nm) paste('-', nm, ' ', options[nm], sep='')), collapse=' '));  
  
  # add gmt files
  gmts<-sapply(yml$gmt, function(x) x[1]);
  gmts<-gmts[file.exists(gmts)];
  gmts<-sapply(strsplit(gmts, '/'), function(x) x[length(x)]);
  gmts[1:length(gmts)]<-paste(prefix, gmts, sep='/'); 
  if (length(gmts) == 0) warning('No available .gmt files') else {
    cmmd<-sapply(names(gmts), function(nm) paste(cmmd, '-gmx', gmts[nm], '-rpt_label', nm));
  }
  
  # Output folder
  cmmd<-paste(cmmd, '-out', prefix);

  # Run GSEA
  if (execute) {
    n<-round(yml$thread[1]);
    if (n>1) {
      run<- parallel::mclapply(cmmd, system, mc.cores = n); 
    } else {
      run<-sapply(cmmd, system);
    }
    
    #file.rename(fn1, paste(out, fn1, sep='/'));
    #fn.out<-paste(prefix, '_', sub('.gmt', '', gmts), '.Gsea', sep='');
    #fn.all<-dir();
    #fn.out<-sapply(fn.out, function(f) fn.all[grep(paste('^', f, sep='')[1], fn.all)]);
    #file.rename(fn.out, paste(out, fn.out, sep='/')); 
    
    SummarizeGSEA(yml$groups[[1]][1], yml$groups[[2]][1], path=prefix, GSEACollection = FALSE);
    file.rename(prefix, paste(out, prefix, sep='/')); 
  }
  

  # save shell command lines
  writeLines(paste(cmmd, '\n\n', sep=''), paste(out, prefix, 'RunGSEA.sh', sep='/'));
  writeLines(yaml::as.yaml(yml), paste(out, prefix, paste(prefix, '.yml', sep=''), sep='/')); 

  cmmd;
  
}