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
# To load this example, download it to your working folder and call:
inputs<-readRDS('./inputs.rds');

# Create report by calling, where _inputs_ is the variable containing all the input data
CreateDeReport(inputs); 

# This version currently requires a collection of gene sets, as in this example: https://github.com/zhezhangsh/DEGandMore/blob/master/examples/DeReport/default_set_human_5-1000.rds?raw=true
# To load this example, download it to your working folder and call:
geneset<-readRDS('./default_set_human_5-1000.rds');
```

##### The inputs is a list with the following elements:
- **anno:** Gene annotation information as a data.frame.
  - the row names are unique Entrez gene IDs
  - the first column must be official gene symbol
	- the last column must be full gene name

 - **expr:** A data matrix includes processed gene expression data.
  - each row corresponds to a gene and each column corresponds to a sample
	- rows must exactly match the rows of **anno**, with the same IDs and the same number of rows
  - it is assumed that the data have been normalized and in linear scale (usually by log2-transformation)

 - **indexes:** Column indexes or column names of samples in **expr**, corresponding to the 2 sample groups. Each group must have more than 1 sample.
  
 - **names:** Names of two sample groups to be compared. By default, the first group is the control or reference.

 - **genome:** Name of the reference genome, in the form of "human", "hsa", or "hg38".

 - **paired:** Boolean indicates if the samples in the 2 groups are paired. Default is **FALSE**. If **TRUE**, the samples will be paired by their order. 
  
 - **penalty:** Whether to give penalty to genes with high sample-sample variance within groups. No penalty if **0**; highest possible penalty if **1**.

 - **homolog:** A list of species to species mapping, that match each gene in **anno** to a gene in another species.

 - **deg:** Options of statistic test to identifiy differentially expressed genes (DEGs).
 
 - **geneset:** Location to file with all the gene sets to be tested for over-representative analsysis of DEGs. 
 
 - **output:** Location to all output files

##### The gene set collection is a list with 2 elements:

 - **meta:** metadata of the genesets. each row is a gene set. include the following columns
   - **Source:** where the gene sets were collected, such as BioSystems, MSigDB, etc.
   - **Collection:** a specific category of the gene sets, such as GO_MF, KEGG_pathway, etc.
   - **Name:** common name of the gene set
   - **Size:** total number of genes in the geneset
   - **URL:** link to the source of the gene set
 - **list:** list of gene sets named by the row names of **meta**, each element corresponds to a row in **meta** and include unique Entrez IDs of genes belonging to the gene set

## Run gene clustering analysis using the _ClReport_ template


