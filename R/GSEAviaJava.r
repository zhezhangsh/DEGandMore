# Run GSEA by creating a Java shell command line
GSEAviaJava<-function(fn.yaml, execute=TRUE) {
  # fn.yaml   Config file specifying inputs
  # execute   Whether to make the GSEA run or just save the command lines
  
  library(yaml);
  yaml<-yaml.load_file(fn.yaml);
  prefix<-yaml$name;
  
  # build shell command line
  cmmd<-paste('java -Xmx5G -cp', yaml$jar);
  
  # output folder
  out<-yaml$out; # folder of run outputs
  if (is.null(out)) out<-'./';
  if (!file.exists(out)) dir.create(out, recursive=TRUE);
   
  # Copy input files to output folder
  fn0<-unlist(yaml[c('input', 'class', 'chip', 'gmt')]);
  fn1<-sapply(strsplit(fn0, '/'), function(fn) fn[length(fn)]);
  file.copy(fn0, paste(out, fn1, sep='/'));
  
  # specify inputs
  if (yaml$preranked) {
    cmmd<-paste(cmmd, 'xtools.gsea.GseaPreranked', '-rnk', fn1['input']);
  } else {
    cmmd<-paste(cmmd, 'xtools.gsea.Gsea', '-res', fn1['input'], '-cls', fn1['class']);
  }
  
  # add chip annotation
  if (!is.null(yaml$chip) & file.exists(yaml$chip)) {
    cmmd<-paste(cmmd, '-chip', fn1['chip']);
  } 
  
  # extra options
  options<-yaml$options;
  cmmd<-paste(cmmd, paste(sapply(names(options), function(nm) paste('-', nm, ' ', options[nm], sep='')), collapse=' '));  
  
  # add gmt files
  gmts<-sapply(yaml$gmt, function(x) x[1]);
  gmts<-gmts[file.exists(gmts)];
  gmts<-sapply(strsplit(gmts, '/'), function(x) x[length(x)]);
  if (length(gmts) == 0) warning('No available .gmt files') else {
    cmmd<-sapply(names(gmts), function(nm) paste(cmmd, '-gmx', gmts[nm], '-rpt_label', paste(prefix, nm, sep='_')));
  }
  
  # Output folder
  cmmd<-paste(cmmd, '-out', './');
  cur.wd<-getwd();
  setwd(out);
  
  # Run GSEA
  if (execute) {
    n<-round(yaml$thread[1]);
    if (n>1) {
      library(snow);
      cl<-makeCluster(n, type='SOCK');
      run<-clusterApplyLB(cl, cmmd, system);
      stopCluster(cl);
    } else {
      run<-sapply(cmmd, system);
    }
    gnm<-paste('Higher_in_', c(yaml$groups[[1]], yaml$groups[[2]]), sep='');
    SummarizeGSEA(gnm[1], gnm[2]);
  }
 
  # save shell command lines
  writeLines(paste(cmmd, '\n\n', sep=''), './RunGSEA.sh');
    
  setwd(cur.wd);
  
  cmmd;
}