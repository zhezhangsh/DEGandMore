# Draw a heatmap with default setting
# 1. in the input data matrix, each row is a variable and each column is a sample
# 2. the input data matrix will be scaled for each row
# 3. columns and rows will be labels by column and row names of the data matrix
# 4. both rows and columns will be clustered using hierarchical clustering 
DegHeatmap<-function(d, col=NA, file.name=NA, cluster.row=TRUE, cluster.col=cluster.row, col.hi='red', col.lo='yellow', plot.new=TRUE) {
library(gplots);

d<-t(scale(t(d)));

if (cluster.row) Rowv<-TRUE else Rowv<-NA;
if (cluster.col) Colv<-TRUE else Colv<-NA;

if (!is.na(file.name) & file.name!='') pdf(paste(file.name, '.pdf', sep=''), w=10, h=10) else if (plot.new) quartz(w=10, h=10);

# set label size
w0<-0.8/max(strwidth(colnames(d), unit='in'));
cexCol<-min(1.5, w0);
w<-0.8/max(strwidth(rownames(d), unit='in'));
h<-(par()$fin[2]-2.5)/nrow(d)/max(strheight(rownames(d), unit='in'));
cexRow<-min(1.5, min(w, h));


# Sample class not specified, plot the heatmap without column side bar; otherwise, color-code the columns
if (identical(col, NA)) heatmap(d, Colv=Colv, Rowv=Rowv, scale='none', cexCol=cexCol, cexRow=cexRow) 
else heatmap(d, scale='none', cexCol=cexCol, Colv=Colv, Rowv=Rowv, cexRow=cexRow, ColSideColors=col[1:ncol(d)]);

if (!is.na(file.name) & file.name!='') dev.off();
}