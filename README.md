# DEGandMore
Analyses related to differential gene expression and functional annotation 

#### Install and load library
```
library(devtools);
install_github("zhezhangsh/DEGandMore");
library(DEGandMore);
```

#### Run GSEA on MSigDB collections and summarize outputs

- Make sure both of Rscript and java are in your PATH
- This example run will save all input and output files into the current directory. Create a new directory and set it as the current directory if you want to
- Copy these input files into your current directory. These files are for the example run only. Prepare your own input file to make a customized run

```
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/gsea/GSEA_example.r               
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/gsea/GSEA_example.yaml              
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/gsea/GSEA_example.cls               
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/gsea/GSEA_example.gct  
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/gsea/gsea2-2.2.0.jar          
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/gsea/h.all.v5.0.symbols.gmt
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/gsea/c2.cp.kegg.v5.0.symbols.gmt 
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/gsea/c2.cp.reactome.v5.0.symbols.gmt
wget https://github.com/zhezhangsh/DEGandMore/raw/master/examples/gsea/GENE_SYMBOL.chip    
```

- Run the GSEA_example script on shell

```
Rscript ./GSEA_example.r ./GSEA_example.yaml

# Alternatively, run the script as below if want to make the actual GSEA run later. 
# The shell script will be save in ./RunGSEA.sh
Rscript ./GSEA_example.r ./GSEA_example.yaml norun

```

- Or, run the example within your R console or RStudio
```
# Only the first time 
library(devtools);
install_github("zhezhangsh/DEGandMore");

# For each run
library(DEGandMore);
GSEAviaJava('GSEA_example.yaml');

```

- Open the index.html file to view results


#### Run pairwise comparison using the _DeReport_ template

```
# Install and load the DEGandMore package
devtools::install_github("zhezhangsh/DEGandMore");
library(DEGandMore);

# Prepare the inputs following this example: https://github.com/zhezhangsh/DEGandMore/blob/master/examples/DeReport/inputs.rds?raw=true

# Create report by calling, where _inputs_ is the variable containing all the input data
CreateDeReport(inputs); 

# This version currently requires a collection of gene sets, as in this example: https://github.com/zhezhangsh/DEGandMore/blob/master/examples/DeReport/default_set_human_5-1000.rds?raw=true

```
