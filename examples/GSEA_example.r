library(devtools);
install_github("zhezhangsh/DEGandMore");
library(DEGandMore);

args<-commandArgs(TRUE);

execute<-TRUE;
if (length(args)>1) if (identical(args[[2]], 'norun')) execute<-FALSE;

cmmd<-GSEAviaJava(args[[1]], execute);
