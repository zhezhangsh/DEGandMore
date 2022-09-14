# Write one or more data matrixes to an Excel file
WriteExcel<-function(values, fileName=paste(Sys.Date(), 'R2Excel', sep='_'), sheetNames=names(values), type='xlsx', settings=rep(list(), length(values)), verbose=TRUE) {
# values		a list of data.frame or matrix
# fileName		the prefix of output file
# sheetNames	the names of individual worksheets
# settings		a list of extra formatting parameters used by WriteExcelSheet, each element of the list corresponds to a worksheet
#               example: rep(list(row.names=FALSE, format=data.frame(ind=3:5, fmt=c('0', '0.00', '0.0E+0'))), 5)
options(java.parameters = "-Xmx102400m");
  
library(xlsx);

#.jinit(classpath="myClasses.jar", parameters="-Xmx512m");

# errorproof for sheet names
if (length(sheetNames)<length(values)) sheetNames[(length(sheetNames)+1):length(values)]<-(length(sheetNames)+1):length(values);
# errorproofing for sheet format setting
if (length(settings)<length(values)) settings[(length(settings)+1):length(values)]<-lapply((length(settings)+1):length(values), function(i) list());
# errorproof for file extension
ext<-paste('\\.', type, '$', sep='');
if (regexpr(ext, fileName, ignore.case=T)==-1) fileName<-paste(fileName, '.', type, sep='');

wb<-createWorkbook(type=type);

#parameters<-lapply(1:length(values), function(i, wb, v, nm, st) append(list(wb=wb, values=v[[i]], name=nm[i]), st[[i]]), wb=wb, v=values, nm=as.character(sheetNames), st=settings);

for (i in 1:length(values)) {
	do.call(WriteExcelSheet, append(list(wb=wb, values=values[[i]], name=as.character(sheetNames[i])), settings[[i]])); 
	if (verbose) print(paste('Create worksheet:', sheetNames[i]));
	try(getSheets(wb), silent=TRUE); # this is a dummy step due to a possible Java bug, as every time the Workbook was modified, the access to this object returns an error
}


saveWorkbook(wb, file=fileName);
} # end of function WriteExcel

####################################################################################################
# Create and format a worksheet in a given Excel Workbook
# Takes some time if the data matrix is large (about 1 min for a 10000X8 matrix)
WriteExcelSheet<-function(wb, values, name, col.names=TRUE, row.names=TRUE, id.name='ID',
    zoom=125, freezeRow=1, freezeCol=0, fontName='Arial', format=data.frame(ind=0, format='')) {
# wb		the Workbook of the Excel file
# values	a data.frame or matrix to write
# name		name of the worksheet
# col.names	whether to write column names
# row.names	whether to write row names
# zoom		adjust the zooming size
# freezeRow where to put the horizontal "freeze pane" bar, no freezing if 0
# freezeCol where to put the vertical "freeze pane" bar, no freezing if 0
# fontName	name of the character font
# format	cell format, such as '0.0000', '0.0E+00' and '0.00%'. A data frame with 2 fields, the first is the column index and the second is the format name

library(xlsx); 

if(row.names) numCol=ncol(values)+1 else numCol=ncol(values);
if(col.names) numRow=nrow(values)+1 else numRow=nrow(values);
colStart<-numCol-ncol(values)+1;
rowStart<-numRow-nrow(values)+1;
header<-colnames(values);

font<-Font(wb, name=fontName);
format<-format[format[,1]>0&format[,1]<=numCol, ];
fmt<-rep('', ncol(values));
fmt[format[,1]]<-as.vector(format[,2]);
styles<-lapply(fmt, function(fmt, wb, font) CellStyle(wb, font=font, dataFormat=DataFormat(fmt)), wb=wb, font=font);

sh<-createSheet(wb, sheetName=name);
cells<-createCell(createRow(sh, rowIndex=1:numRow), colIndex=1:numCol);

    assign<-function(c, v, s=CellStyle(wb)) {
        sapply(1:length(c), function(i, c, v, s) {setCellStyle(c[[i]], s); setCellValue(c[[i]], v[i]); c[[i]]}, c=c, v=as.vector(v), s=s)
    };

# write column names
if (col.names) {
if (row.names) header<-c(id.name, header);
#cells[1,]<-assign(cells[1,], header, CellStyle(wb, font=Font(wb, name=fontName, isBold=TRUE)));
assign(cells[1,], header, CellStyle(wb, font=Font(wb, name=fontName, isBold=TRUE)));
}
# write row names
#if (row.names) cells[rowStart:numRow,1]<-assign(cells[rowStart:numRow,1], rownames(values), CellStyle(wb, font=font));
assign(cells[rowStart:numRow,1], rownames(values), CellStyle(wb, font=font));

# assign values to cells
#cells[rowStart:numRow, colStart:numCol]<-
#    sapply(1:ncol(values), function(i, c, v, s) assign(c[, i], v[, i], s[[i]]), c=cells[rowStart:numRow, colStart:numCol], v=values, s=styles);
sapply(1:ncol(values), function(i, c, v, s) if (class(c)[1]=='matrix') assign(c[, i], v[, i], s[[i]]) else assign(c[i], v[, i], s[[i]]),
    c=cells[rowStart:numRow, colStart:numCol], v=values, s=styles);

setZoom(sh, numerator=zoom, denominator=100);
if (freezeRow>0|freezeCol>0) createFreezePane(sh, rowSplit=freezeRow+1, colSplit=freezeCol+1);
autoSizeColumn(sh, 1);
}


