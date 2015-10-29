# Generate a GSEA report in HTML format, with hyperlinks to subfolders
SummarizeGSEA<-function(name1, name2, fn.index='index.html', path='.', GSEACollection=TRUE) {
    # name1, name2          Column names
    # title                 Title of the HTML file
    # fn                    Name of the HTML file
    # path                  Working directory, where all the GSEA results locate
    # GSEACollection        Whether the GSEA gene set collections themseleves were used for analysis
    
    wd0<-getwd();
    setwd(path);
  
    # Top level collections
    collections<-c(
       HALLMARK = "Collection 0: Hallmark gene sets",
       C1 = "Collection 1: Positional gene sets",
       C2 = "Collection 2: Curated gene sets",
       C3 = "Collection 3: Motif gene sets",
       C4 = "Collection 4: Computational gene sets",
       C5 = "Collection 5: GO gene sets",
       C6 = "Collection 6: Oncogenic signatures gene sets",
       C7 = "Collection 7: Immunologic signatures gene sets"
    );
    
    # Second level collections
    collections2<-c(
       "h.all" = "Collection 0: Hallmark gene sets",
       "CGP" = "Collection 2: Chemical and genetic perturbations",
       "CP" = "Collection 2: All canonical pathways",
       "BIOCARTA" = "Collection 2: BioCarta pathways",
       "KEGG" = "Collection 2: KEGG pathways",
       "REACTOME" = "Collection 2: Reactome pathways",
       "MIR" = "Collection 3: MicroRNA targets",
       "TFT" = "Collection 3: Transcription factor targets",
       "CGN" = "Collection 4: Cancer gene neighborhoods",
       "CM" = "Collection 4: Cancer modules",
       "BP" = "Collection 5: GO biological processes",
       "CC" = "Collection 5: GO cellular components",
       "MF" = "Collection 5: GO molecular functions"
    );
 
    f<-dir();
    f<-f[file.info(f)[,'isdir']];
    #f<-f[grep('.Gsea.', f)];
    
    # Original summarization file
    pos<-sapply(f, function(f) {
      fn<-dir(f); 
      fn[grepl('^gsea_report_for_', fn) & grepl('.html$', fn) & (grepl('_na_pos_', fn) | grepl(paste('_', name1, '_', sep=''), fn))][1];
    })
    neg<-sapply(f, function(f) {
      fn<-dir(f); 
      fn[grepl('^gsea_report_for_', fn) & grepl('.html$', fn) & (grepl('_na_neg_', fn) | grepl(paste('_', name2, '_', sep=''), fn))][1];
    })

    flg<-!is.na(pos) | !is.na(neg);
    pos<-pos[flg];
    neg<-neg[flg];
    f<-f[flg];
    
    # create the index.html file
    if (length(pos)>0 & length(pos)==length(neg)) {
      nm<-sapply(strsplit(f, '.Gsea.'), function(x) x[1]);    
      file.rename(f, nm); 
      pos<-paste(nm, pos, sep='/');
      neg<-paste(nm, neg, sep='/');
      if (GSEACollection) {
        for (i in 1:length(collections2)) nm[grepl(names(collections2)[i], nm, ignore=TRUE)][1]<-collections2[i];       
        for (i in 1:length(collections)) nm[grepl(toupper(names(collections)[i]), toupper(nm))][1]<-collections[i];
      }
          
      urls<-cbind(pos, neg);
      rownames(urls)<-nm; 
      urls<-urls[order(rownames(urls)), , drop=FALSE]; 
      urls[, 1]<-awsomics::AddHref('List', urls[, 1]);
      urls[, 2]<-awsomics::AddHref('List', urls[, 2]);
      urls<-cbind(rownames(urls), urls);
      colnames(urls)<-c('Gene set collection', paste(name1, name2, sep=' > '), paste(name2, name1, sep=' > '));
      awsomics::CreateDatatable(urls, fn.index, rownames = FALSE, caption = paste('GSEA Result Index:', name1, 'vs.', name2));
       
      # Files with full table of results
      fn1<-sub('.html$', '.xls', c(pos, neg));
      nm1<-rep(sub('Collection ', '', nm), 2);
      tbls<-lapply(fn1, function(f) read.csv2(f, sep='\t', row=1));
      coll<-rep(nm1, as.vector(sapply(tbls, nrow)));
      tbl<-do.call('rbind', tbls);
      tbl<-data.frame(Collection=coll, Gene_set=tbl[[1]], Size=as.vector(tbl[[3]]), NES=as.vector(tbl[[5]]), PValue=as.vector(tbl[[6]]), FDR=as.vector(tbl[[7]]), stringsAsFactors=FALSE);
      for (i in 3:6) tbl[, i]<-as.numeric(tbl[, i]);
      tbl<-tbl[order(tbl$NES), ];
      di<-rep('No', nrow(tbl));
      di[tbl$NES<0]<-'Yes';
      tbl<-data.frame(tbl, di, stringsAsFactors=FALSE);
      names(tbl)[ncol(tbl)]<-paste(name2, name1, sep='>');
      saveRDS(tbl, file='full_table.rds');
      tbl$Gene_set<-awsomics::AddHref(tbl$Gene_set, paste("http://www.broadinstitute.org/gsea/msigdb/cards", tbl$Gene_set, sep='/'));
      awsomics::CreateDatatable(tbl, "full_table.html", rownames = FALSE, caption = paste('GSEA Result Table:', name1, 'vs.', name2));
    }
    
    setwd(wd0);
}

