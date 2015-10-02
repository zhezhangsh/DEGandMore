
# Export functions in NAMESPACE file
path<-'./R';
fn<-paste(path, dir(path), sep='/');
fn<-fn[grep('.r$', fn, ignore.case=TRUE)];
if (length(fn)>0) {
	fnc<-lapply(fn, function(fn) {
		ln<-scan(fn, what='', flush=TRUE, sep='\n');
    ln<-ln[!grepl('^ ', ln)];
		ln<-gsub(' ', '', ln);
		ln<-ln[grep("<-function\\(", ln)];
		sapply(strsplit(ln, "<-function\\("), function(x) x[1]);
	});
	fnc<-sort(as.vector(unlist(fnc)));
  ln<-sort(paste('export("', fnc, '");', sep=''));
  writeLines(c('', ln, ''), './NAMESPACE');
}

