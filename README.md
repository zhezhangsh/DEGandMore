# DEGandMore
Analyses related to differential gene expression and functional annotation 

#### Install and load library
```
library(devtools);
install_github("zhezhangsh/DEGandMore");
library(DEGandMore);
```

#### Make an example run of GSEA via command line

- Make sure both of Rscript and java are in your PATH
- This example run will save all input and output files into the current directory. Create a new directory and set it as the current directory if you want to
- Copy these input files into your current directory. These files are for the example run only. Prepare your own input file to make a customized run
```
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/GENE_SYMBOL.chip
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/
```
- Run the GSEA_example script
```
Rscript ./
```
- Alternatively, run below if want to make the actual GSEA run later. The shell script will be save in ./RunGSEA.sh
```

```
