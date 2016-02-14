# Given a table of gene list, write a detailed table with hyperlinks to an HTML
GeneList2Datatable<-function(d, fn='gene_list', col.entrez=0, col.symbol=1, genome=NA, title='Gene_List', extra.url=list()) {
    # d     A data.frame of gene list, with entrez ID as rowname by default, and gene symbol in the first column. The other columns could be more annotation or statistics related to each gene
    # col.entrez, col.symbol  The column number of key ID types. -1: not available; 0: row name
    # extra.url     Extra column(s) to be attached to hyperlinks; elements of list should be named after column names and each element has the value of hyperlinks
  
    d<-as.data.frame(d);

    # HTML file name
    fn<-sub('.html$', '', fn, ignore.case=TRUE);
    fn<-sub('.htm$', '', fn, ignore.case=TRUE);
    fn<-paste(fn, '.html', sep='');

    if (nrow(d)==0) {        
      CreateDatatable(d, fn, rownames=FALSE, caption=paste(title, '(Empty table)'));
    } else {
    
    if (is.na(col.entrez)) col.entrez<--1;
    if (is.na(col.symbol)) col.symbol<--1;
    if (col.entrez==0 | col.symbol==0) {
        d<-cbind(ID=rownames(d), d);
        col.entrez<-col.entrez+1;
        col.symbol<-col.symbol+1;
    }
    
    # build UCSC url
    cn<-col.symbol;
    if (cn<0 | is.na(cn) | cn>ncol(d)) ucsc<-NA else {
        ucsc<-paste('http://genome.ucsc.edu/cgi-bin/hgTracks?position=', as.vector(d[, col.symbol]), sep=''); 
        if (!is.na(genome)) {
            gnm<-awsomics::GetGenomeAlias(genome, 'ucsc');
            if (!is.na(gnm)) ucsc<-paste(ucsc, '&db=', gnm, sep='');
        }
    }
    
    # build KEGG url
    cn<-col.entrez;
    if (cn<0 | is.na(cn) | cn>ncol(d)) kegg<-NA else {
        #if (cn==0) id<-rownames(d) else id<-as.vector(d[, cn]);
        kegg<-paste('http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&keywords=%3A', as.vector(d[, col.entrez]), sep='');
        if (!is.na(genome)) {
            gnm<-awsomics::GetGenomeAlias(genome, 'code');
            if (!is.na(gnm)) {
                if (gnm=='cel') {
                    if (col.symbol==-1) kegg<-NA # KEGG use gene symbol for worms
                    else kegg<-paste('http://www.genome.jp/dbget-bin/www_bget?', gnm, ':', as.vector(d[, col.symbol]), sep='') 
                } else kegg<-paste('http://www.genome.jp/dbget-bin/www_bget?', gnm, ':', as.vector(d[, col.entrez]), sep='');
            }
        }
    }
    
    # build ENTREZ url
    cn<-col.entrez;
    if (cn<0 | is.na(cn) | cn>ncol(d)) entrez<-NA else entrez<-awsomics::UrlEntrezGene(d[, col.entrez]);

    if (!identical(NA, kegg)) { # add column to KEGG pathway
        d$Pathway<-rep('KEGG', nrow(d));
        d<-d[, c(1, ncol(d), 2:(ncol(d)-1))];
    }
    if (!identical(NA, ucsc)) { # add column to Genome Browser
        d$Browser<-rep('UCSC', nrow(d));
        d<-d[, c(1, ncol(d), 2:(ncol(d)-1))];
    }
    
    # extra column(s) to add hyperlinks
    extra.url<-extra.url[names(extra.url) %in% colnames(d)];
    if (length(extra.url)>0) {
        cnm<-colnames(d);
        ex<-sapply(extra.url, function(x) x[rownames(d)]);
        for (i in 1:ncol(ex)) { 
            nm<-colnames(ex)[i];
            d<-transform(d, temp=paste('<a href = ', shQuote(ex[,i]), '>', as.vector(d[,nm]), '</a>', sep=''));
            d[, nm]<-as.vector(d[, 'temp']);
            d<-d[, colnames(d)!='temp'];
            d[ex[,nm]=='' | is.na(ex[,nm]),nm]<-'';
        }
        colnames(d)<-cnm;
    }
    
    if (!identical(NA, entrez)) colnames(d)[col.entrez]<-'Entrez_ID';
    nm<-colnames(d);
    if (!identical(NA, entrez)) d$Entrez_ID<-awsomics::AddHref(d$Entrez_ID, entrez);
    if (!identical(NA, ucsc)) d$Browser<-awsomics::AddHref(d$Browser, ucsc);
    if (!identical(NA, kegg)) d$Pathway<-awsomics::AddHref(d$Pathway, kegg);

    d$Entrez_ID[is.na(entrez) | entrez=='']<-'';
    d$Browser[is.na(ucsc) | ucsc=='']<-'';
    d$Pathway[is.na(kegg) | kegg=='']<-'';
    colnames(d)<-nm; # prevent column names to be changed by the calls

    title<-gsub(' ', '_', title);
    awsomics::CreateDatatable(d, fn, rownames=FALSE, caption=title);
   } # end of nrow(d)>0

d;
}