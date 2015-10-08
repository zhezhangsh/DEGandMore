# Run a GSEA analysis of KEGG pathways on pairwise comparison results 
# Pathway genes: KEGG database
# Differential expression: RankProd
# Gene set analysis: GAGE
# Result illustration: Pathview

GseaKegg<-function(e, g1.ind, g2.ind, groups=c('A', 'B'), paired=FALSE, genome='hsa', path='.', at.least.plot=5, at.most.plot=Inf, 
                      ranking.penalty=0, num.thread=1, path.xml='.') {
    # e                 Data matrix of input data, row names must be Entrez gene ID
    # g1.ind, g2.ind    Column indexes or names of 2 groups of samples in e
    # groups            Names of the two compared groups
    # paired            Whether it's paired comparison
    # genome            Genome name, assembly, or code, such as 'worm', 'mm10', and 'human'
    # path              Path to output files
    # at.least.plo      Minimum number of pathways to plot
    # ranking.penalty   The penalty on high-variance genes during gene ranking. Value between 0 to 1. No penalty if 0. Larger value means more penalty on high variance genes
    # num.thread        Number of threads to use for running pathview function in parallel

    library(org.Hs.eg.db);
    library(gage);
    library(pathview);
    library(RankProd)
    library(DEGandMore);    
    
    if (at.least.plot<1) at.least.plot<-1; # at least one pathway to plot
    
    g1.name<-groups[1];
    g2.name<-groups[2];
       
    g1.ind<-colnames(e[, g1.ind, drop=FALSE]);
    g2.ind<-colnames(e[, g2.ind, drop=FALSE]);
    
    genome<-GetGenomeAlias(genome, 'code'); # get genome code for querying KEGG
    
    if (is.na(genome)) {
        cat('Unknown genome', genome, '\n');
        kegg<-NA;
    } else {
        gs<-gage::kegg.gsets(genome, id.type='entrez')[[1]];
        
        # Mapping worm IDs because KEGG doesn't use ENTREZ ID for worm genes
        #if (genome == 'cel') {
        #    gs.cel<-gs;
        #    mp<-lapply(gs, function(gs) id2eg(sub('^CELE_', '', gs), 'ENSEMBL', pkg.name=SC.GetGenomeAlias('cel', 'package')));
        #    gs<-lapply(mp, function(mp) mp[,'EntrezGene'], '; '))));
        #}
        
        # Run gene set analysis with GAGE
        ############################################################################################################
        ############################################################################################################
        if (paired) pr<-'paired' else pr<-'unpaired';
        #g1<-colnames(e[, g1.ind]);
        #g1<-(1:ncol(e))[which(colnames(e) %in% g1)]
        #g2<-colnames(e[, g2.ind]);
        #g2<-(1:ncol(e))[which(colnames(e) %in% g2)]

        gg<-gage::gage(e[, c(g1.ind, g2.ind)], gs, ref=1:length(g1.ind), samp=(length(g1.ind)+1):(length(g1.ind)+length(g2.ind)), same.dir=TRUE, compare=pr);
        stat.gg<-lapply(gg[1:2], function(x) x[, c(5, 2, 3, 4)]);
        sig.gg<-lapply(stat.gg, function(x) {
            ind<-which(!is.na(x[,3]) & (x[,3]<=0.05 | x[,4]<=0.25));
            if (length(ind)<at.least.plot) x[1:min(at.least.plot, nrow(x)), ] else x[ind, , drop=FALSE];
        });
        sig.gg<-lapply(sig.gg, function(x) x[1:min(at.most.plot, nrow(x)), , drop=FALSE]);
        Gage<-list(significant=sig.gg, stat=stat.gg, full=gg);
        ############################################################################################################

        # Run differential expression with Rank Product 
        ############################################################################################################
        ############################################################################################################
        if(ranking.penalty>0) {
            if (penalty>1) penalty<-1;
            means<-cbind(rowMeans(e[, g1.ind]), rowMeans(e[, g2.ind]));
            l2r<-means[,2]-means[,1];
            fc<-exp(l2r*log(2));
            colnames(means)<-paste('Mean_', c(g1.name, g2.name), sep='');
            if (paired)
                sd<-apply(e[, g2.ind]-e[, g1.ind], 1, sd)
            else {
                v1<-apply(e[, g1.ind], 1, var);
                v2<-apply(e[, g2.ind], 1, var);
                df1<-length(g1.ind)-1;
                df2<-length(g2.ind)-1;
                sd<-sqrt((df1*v1+df2*v2)/(df1+df2));
                sd[sd==0]<-min(sd[sd>0])/2;
            }
            pnl<-quantile(sd, probs=seq(0, 1, 0.01))[100*round(1-ranking.penalty,2)+1];
            temp<-e; # save the original
            e[, g2.ind]<-apply(e[, g2.ind], 2, function(d) means[,1]+(d-means[,1])/(pnl+sd));
        }
        
        if (paired) {
            d<-e[, g1.ind]-e[, g2.ind];
            cl<-rep(1, ncol(d));
        } else {
            d<-e[, c(g1.ind, g2.ind)];
            cl<-rep(0:1, c(length(g1.ind), length(g2.ind)));
        }
        rp<-RankProd::RP(d, cl);
        rp<-SummarizeRP(rp, groups[1], groups[2], rownames(d), FALSE, TRUE, FALSE, FALSE);
        rp0<-log2(1/rp[[1]][,2])+log2(rp[[2]][,2]);
        if (ranking.penalty>0) {
            rp[[3]][,4]<-rp[[1]][,5]<-rp[[2]][,5]<-l2r;
            rp[[3]][,5]<-rp[[1]][,6]<-rp[[2]][,6]<-fc;
            e<-temp;
        }
        stat.rp<-cbind(rp[[3]], RP=rp0);
        RankProd<-list(stat=stat.rp, full=rp);
        ############################################################################################################
        
        # Plot pathways via pathview
        ############################################################################################################
        ############################################################################################################
        ids<-lapply(sig.gg, function(x) sapply(strsplit(rownames(x), ' '), function(x) x[1]));
        ids<-sort(unique(unlist(ids)));
        v<-sort(stat.rp[,4])[rank(stat.rp[,8])];
        
        if (file.exists(path.xml)) {
          fn.xml<-dir(path.xml);
          fn.xml<-fn.xml[grep(genome, fn.xml)];
          fn.xml<-fn.xml[grep('.xml$', fn.xml)];
          nm<-sub('.xml$', '', fn.xml);
          fn.xml<-paste(path.xml, fn.xml, sep='/');
          names(fn.xml)<-nm;
          fn.xml<-fn.xml[ids];
          names(fn.xml)<-ids;
          fn.xml[is.na(fn.xml)]<-'';
        } else {
          fn.xml<-rep('', length(ids));
          names(fn.xml)<-ids;
        }
        
        # file folder skeleton
        if (!file.exists(path)) dir.create(path);
        pth.rs<-paste(path, '/Files', sep='');        
        if (!file.exists(pth.rs)) dir.create(pth.rs);
        figf<-paste(pth.rs, '/figures', sep='');
        if (!file.exists(figf)) dir.create(figf)
        subf<-paste(path, '/', c("Higher_in_", "Lower_in_"), g2.name, sep='');
        sapply(subf, function(subf) if (!file.exists(subf)) dir.create(subf));
        sb<-paste(figf, 'colored_plots', sep='/');
        if(!file.exists(sb)) dir.create(sb, recursive = TRUE);
        sb.genes<-paste(subf, '/genes', sep='');
        sapply(sb.genes, function(f) if (!file.exists(f)) dir.create(f));
        sb.ht<-paste(figf, '/heatmaps', sep='');
        sapply(sb.ht, function(f) if (!file.exists(f)) dir.create(f));
        sb.bar<-paste(subf, '/bars', sep='');
        sapply(sb.bar, function(f) if (!file.exists(f)) dir.create(f));
        
        wd<-getwd();
        setwd(sb);
        #xmls<-paste(kegg.dir, '/', id, '.xml', sep='');
        #has.local<-file.exists(xmls);
        if (num.thread > 1) {
          pathviewCl<-function(ids, gene.data=v, species=genome, xml=fn.xml, dir=path.xml, low=list(gene = "blue", cpd = "green")) {
            library(pathview);
            lapply(ids, function(id) pathview(gene.data, pathway.id=id, xml.file=xml[id], kegg.dir=dir, species=species, low = low));
          }
          
          ids<-split(ids, cut(1:length(ids), num.thread, labels=FALSE));
          library(snow);
          cl<-makeCluster(num.thread, type='SOCK');
          pv<-clusterApplyLB(cl, ids, pathviewCl, gene.data=v, species=genome, xml=fn.xml, dir=path.xml); 
          pv<-unlist(pv, recursive=FALSE);
          ids<-unlist(ids, use.names=FALSE);
        }
        else {
          capture.output(pv<-lapply(ids, function(id) pathview(v, pathway.id=id, species=genome, low = list(gene = "blue", cpd = "green"), xml.file=fn.xml[id], kegg.dir=path.xml)))->x;
        } 
        setwd(wd);
        
        names(pv)<-ids;
        ############################################################################################################
        
        # Writing results to files
        ############################################################################################################
        ############################################################################################################
        ############################################################################################################        
        
        # copy plots to subfolders
        ############################################################################################################        
        ############################################################################################################        
#         lapply(ids, function(id) {
#             fn<-paste(path, paste(id, c('.xml', '.png', '.pathview.png'),sep=''), sep='/'); 
#             file.rename(fn[file.exists(fn)], paste(sb[file.exists(fn)], fn[file.exists(fn)], sep='/'));
#         })
        ############################################################################################################
       
        # Write gene list of each pathway
        ############################################################################################################        
        ############################################################################################################        
        # Get gene symbol
        mp<-eg2id(rownames(e), 'SYMBOL', pkg.name=awsomics::GetGenomeAlias(genome, 'package'));
        symb<-mp[,2];
        names(symb)<-mp[,1];
        symb<-symb[rownames(e)];
        names(symb)<-rownames(e);
        
        # Get pathway tables
        tbls<-lapply(1:2, function(i) {
            stat<-rp[[i]];
            stat<-data.frame(Symbol=symb[rownames(stat)], Rank=stat[,1], rowMeans(e[, g1.ind]), rowMeans(e[, g2.ind]), stat[, 5],Fold_Change=stat[,6],P_Value=stat[,3],FDR=stat[,4]);
            colnames(stat)[3:5]<-c(paste('Mean_', g1.name, sep=''), paste('Mean_', g2.name, sep=''), paste(g2.name, '-', g1.name, sep=''));
            stat<-stat[order(stat[,'Rank']), ];
            tbl<-lapply(gs, function(gs) stat[rownames(stat) %in% gs, , drop=FALSE]);
            names(tbl)<-sapply(strsplit(names(gs), ' '), function(x) x[1]);
            tbl;
        });
        
        gs.nm<-names(gs);
        names(gs.nm)<-sapply(strsplit(names(gs), ' '), function(x) x[1]);

        # writing tables to html files
        lapply(1:2, function(i) {
            t<-tbls[[i]];
            sapply(names(t), function(nm) GeneList2Datatable(awsomics::FormatNumeric(t[[nm]]), fn=paste(sb.genes[i], nm, sep='/'), genome=genome, title=nm));
        })->x;
        ############################################################################################################

        # Plot heatmap of all pathways
        ############################################################################################################
        ############################################################################################################
        rp0<-stat.rp;
        rp0<-rp0[order(rp0[,1]),];
        t<-lapply(names(gs), function(nm) {
            ht<-e[rownames(e) %in% gs[[nm]], , drop=FALSE];
            col<-rep(c('#888888', '#3333FF'), c(length(g1.ind), length(g2.ind)));
            if (nrow(ht)>=2) GseaHeatmap(ht, as.vector(RankProd[[1]][rownames(ht), 1]), g1.ind, g2.ind, title=nm, file.name=paste(sb.ht, '/', strsplit(nm, ' ')[[1]][1], sep=''));
            #if (nrow(ht)>=2)SC.Heatmap(ht, col=col, cluster.row=FALSE, cluster.col=TRUE, file.name=paste(sb.ht, '/', strsplit(nm, ' ')[[1]][1], sep=''));
            ht[, c(g1.ind, g2.ind), drop=FALSE];
        });
        ############################################################################################################        
        
        # Plot bar distributions of all pathways
        ############################################################################################################
        ############################################################################################################   
        rk<-lapply(1:2, function(i) {
            rp0<-rp[[i]];
            if (i==1) rv<-FALSE else rv<-TRUE;
            lapply(names(tbls[[i]]), function(nm) { 
                rp1<-rp0[rownames(rp0) %in% rownames(tbls[[i]][[nm]]), , drop=FALSE];
                ht<-rp1[, 'Rank'];
                names(ht)<-rownames(rp1);
                fn<-paste(sb.bar[i], '/', nm, sep='');
                RunningRank(ranks=ht, n=nrow(rp0), reverse=rv, fn=fn, title=gs.nm[nm]);
            });
        });
        ############################################################################################################        
        

        # output pathway enrichment tables
        ############################################################################################################        
        ############################################################################################################        
        formatted<-lapply(1:2, function(i) {
            gg<-stat.gg[[i]];
            id<-sapply(strsplit(rownames(gg), ' '), function(x) x[1]);
            nm<-substr(rownames(gg), nchar(id)+2, nchar(rownames(gg)));
            #n<-sapply(gs, function(gs) length(intersect(gs, rownames(e))))[rownames(gg)];
            n<-sapply(gs, length)[rownames(gg)];
            tbl<-data.frame(Pathway_ID=id, Num_Genes_Total=n, Num_Genes_Measured=gg[, 'set.size'], Enrichment_Score=gg[, 'stat.mean'], P_Value=gg[, 'p.val'], FDR=gg[, 'q.val'], Genes=rep('List', nrow(gg)), Ranking=rep('Figure', nrow(gg)), Heatmap=rep('Figure', nrow(gg)),  Pathview=rep('Figure', nrow(gg)), Name=nm);
            #tbl<-data.frame(Pathway_ID=id, Num_Genes_Total=n, Num_Genes_Measured=gg[, 'set.size'], Enrichment_Score=gg[, 'stat.mean'], P_Value=gg[, 'p.val'], FDR=gg[, 'q.val'], Ranking=rep('Figure', nrow(gg)), Heatmap=rep('Figure', nrow(gg)),  Pathview=rep('Figure', nrow(gg)), Name=nm);
            tbl.raw<-tbl;
            fmt<-list(Enrichment_Score='0.00', P_Value='0.00E+0', FDR='0.00%');
            
            url.kegg<-paste('http://www.genome.jp/kegg-bin/show_pathway?', id, sep='');
            url.list<-paste('./genes/', id, '.html', sep='');
            url.rank<-paste('./bars/', id, '.pdf', sep='');
            url.heat<-paste('../Files/figures/heatmaps/', id, '.pdf', sep='');
            url.view<-paste('../Files/figures/colored_plots/', id, '.pathview.png', sep='');

            tbl<-transform(tbl, 'Pathway_ID'=paste('<a href = ', shQuote(url.kegg), '>', as.vector(tbl[,'Pathway_ID']), '</a>', sep=''));
            tbl<-transform(tbl, 'Genes'=paste('<a href = ', shQuote(url.list), '>', as.vector(tbl[,'Genes']), '</a>', sep=''));
            tbl<-transform(tbl, 'Ranking'=paste('<a href = ', shQuote(url.rank), '>', as.vector(tbl[,'Ranking']), '</a>', sep=''));
            tbl<-transform(tbl, 'Heatmap'=paste('<a href = ', shQuote(url.heat), '>', as.vector(tbl[,'Heatmap']), '</a>', sep=''));
            tbl<-transform(tbl, 'Pathview'=paste('<a href = ', shQuote(url.view), '>', as.vector(tbl[,'Pathview']), '</a>', sep=''));
            
            tbl[!file.exists(paste(figf, '/heatmaps/', id, '.pdf', sep='')), 'Heatmap']<-'';
            tbl[!file.exists(paste(figf, '/colored_plots/', id, '.pathview.png', sep='')), 'Pathview']<-'';

            if(!file.exists(subf[i])) dir.create(subf[i]);
            awsomics::CreateDatatable(awsomics::FormatNumeric(tbl), fn=paste(subf[i], 'index.html', sep='/'), rownames=FALSE);
            tbl.raw;
        });
        
        # summary table
        n<-sapply(stat.gg, function(x) c(nrow(x[!is.na(x[,3])&x[,3]<=0.05, , drop=FALSE]), c(nrow(x[!is.na(x[,3])&x[,4]<=0.25, , drop=FALSE]))));
        chg<-c(paste(c('Higher', 'Lower'), '_in_', g2.name, sep=''));
        smm<-data.frame(Change=chg, t(n));        
        smm<-transform(smm, 'Change'=paste('<a href = ', shQuote(paste('./', c('Higher_in_', 'Lower_in_'), g2.name, '/index.html', sep='')), '>', as.vector(smm[,'Change']), '</a>', sep=''));
        colnames(smm)[2:3]<-c('p < 0.05', 'FDR < 0.25');
        
        awsomics::CreateDatatable(smm, paste(path, '/index.html', sep=''), rownames = FALSE);
        #print(gvisTable(smm, options = list(allowHTML = TRUE, width='400px', height='200', showRowNumber=TRUE), chartid='Summary_Kegg_Analysis'), file=paste(path, '/index.html', sep='')); 

        kegg<-list(Data=e, Genome=genome, Group_Names=groups, Groups=list(g1.ind, g2.ind), Paired=paired, Gene_sets=gs, Formatted=formatted, GAGE=Gage, RankProd=RankProd, Pathview=pv);
         
        write.csv(kegg$formatted[[1]], paste(pth.rs, '/Pathways_higher_in_', g2.name, '.csv', sep=''));
        write.csv(kegg$formatted[[2]], paste(pth.rs, '/Pathways_lower_in_', g2.name, '.csv', sep=''));
        write.csv(kegg$RankProd[[1]], paste(pth.rs, '/Genes_higher_in_', g2.name, '.csv', sep=''));
        write.csv(kegg$RankProd[[2]], paste(pth.rs, '/Genes_lower_in_', g2.name, '.csv', sep=''));
        saveRDS(kegg, file=paste(pth.rs, '/kegg.rds', sep=''));
    };
       
    kegg;
}


# plot a GSEA style heatmap of genes in a specific pathway
GseaHeatmap<-function(d, rnk, g1.ind, g2.ind, title='', file.name=NA) {
    # d                 The data matrix
    # rnk               Rank of genes to order the matrix, must have the same length the number of rows of "d"
    # g1.ind, g2.ind    Sample index
    # title             Pathway name
    # file.name         File name of output pdf file
    
    d<-d[, c(g1.ind, g2.ind)];
    d<-d[order(rnk), ];
    #dn<-d-rowMeans(d);
    dn<-t(scale(t(d))); # normalize data
    dn[dn < -3]<--3;
    dn[dn > 3]<-3;

    # determin plot size
    W<-0.2*ncol(d) + 0.1*max(nchar(rownames(d))) + 0.5;
    H<-0.2*nrow(d) + 0.1*max(nchar(colnames(d))) + 0.5;
    
    if (!is.na(file.name) & file.name!='') pdf(paste(file.name, '.pdf', sep=''), w=W, h=H+0.5) else if (plot.new) quartz(w=W, h=H+0.5);
    
    par(mar=c(0,0,0,0), omi=c(0.1, 0.1, 0.6, 0.1));
    plot(0, type='n', xlim=c(0, W/0.2-1), ylim=c(0, H/0.2-1), xaxs='i', yaxs='i', axes=FALSE, xlab='', ylab='');
    
    col=rep(c("#CCCCCC", "#FFEE00"), c(length(g1.ind), length(g2.ind)));
    sapply(1:length(col), function(i) rect(i-1, nrow(dn), i, H/0.2-1, border=NA, col=col[i]))->x;

    # blue-pinkogram colors. The first and last are the colors to indicate the class vector (phenotype). This is the 1998-vintage, pre-gene cluster, original pinkogram color map
    col<-c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000");

    image(0:ncol(dn), 0:nrow(dn), t(dn[nrow(dn):1,]), col=col, add=TRUE);
    text(W/0.2-1.5, (1:nrow(dn))-0.5, offset=0, pos=2, label=rownames(dn)[nrow(dn):1]);
    text((1:ncol(dn))-0.5, H/0.2-1.5, pos=2, offset=0, srt=90, label=colnames(dn));
    title(main=title, cex.main=0.8*(W-0.2)/strwidth(title, unit='in'), outer=TRUE);

    abline(v=0:ncol(dn), h=0:nrow(dn), lwd=0.5);
    box();
    dev.off();
}

# Make a GSEA-style plot of running ranks
RunningRank<-function(ranks, n, reverse=FALSE, line.col=if (reverse) 4 else 2, fn=NA, title='') {
  # ranks         Ranks of a subset of elements
  # n             Total number of elements
  # fn            PDF file name to output if not NA
  nm<-names(ranks);
  
  ranks<-ranks[ranks>0&ranks<=n];
  
  if (reverse) ranks<-n+1-ranks;
  
  if (!identical(fn, NA)) pdf(paste(fn, '.pdf', sep=''), w=8, h=4.8);
  
  if (length(ranks)==0) {
    plot.new();
    title(main='No elements', cex.main=2);
  } else {
    ranks<-sort(ranks);
    rt<-(1:length(ranks))/((ranks/n)*length(ranks));
    x<-c(0, ranks, n);
    y<-c(0, log2(rt), 0);
    mx<-max(y);
    mn<-min(y);
    
    if (reverse) {  
      xlim=c(n, 0);
      ylim=c(mn-0.25*(mx-mn)-0.1*(mx-mn), mx);
      smmt<-which(y==min(y));
      smmt<-smmt[length(smmt)];
      smmt.y<-min(y)-0.05*(mx-mn);
      smmt.p<-2;
    } else {
      xlim=c(0, n);
      ylim=c(mn-0.25*(mx-mn), mx+0.1*(mx-mn));
      smmt<-which(y==max(y));
      smmt<-smmt[length(smmt)];
      smmt.y<-max(y)+0.05*(mx-mn);
      smmt.p<-4;
    }
    
    par(mar=c(6,5,3,2));
    plot(x, y, type='l', lwd=2, lty=1, cex=1, col=line.col, xaxs='i', xlim=xlim, ylim=ylim, xlab='Gene ranking index', ylab='Fold enrichment', sub=paste('Number of genes: ', n, ' (total), ', length(x), ' (in gene set)', sep=''), cex.lab=2, main=title, cex.main=2);
    #points(x, y, pch='|', cex=0.5); 
    abline(h=0, lty=3); 
    lb<-par()$usr[3];
    rt<-lb+0.2*(par()$usr[4]-par()$usr[3]);
    #rect(0, lb, n, rt, col='#3333FF88', border=NA); 
    #abline(h=rt);
    
    abline(v=x[smmt], lty=2, col='grey');
    text(x[smmt], smmt.y, pos=smmt.p, label=paste("Peak at gene", nm[smmt]));
    if (reverse) rk<-n+1-ranks else rk<-ranks;
    segments(rk, lb+0.02*(mx-mn), rk, rt, col='darkgrey', lwd=min(1, max(0.3, 100/length(x))));   
  }
  
  if (!identical(fn, NA)) dev.off();
}

