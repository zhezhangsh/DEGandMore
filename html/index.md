---
title: "RNA-seq 2G"
author: "Zhe Zhang"
date: "2016-09-26"
output:
  html_document:
    number_sections: yes
    self_contained: no
    toc: yes
    toc_float:
      collapsed: no
---

<div style="border:black 1px solid; padding: 0.5cm 0.5cm"><font color="darkblue">
**RNA-seq 2G** is a web portal with >20 statistical methods that perform two-group analysis of differential gene expression. It uses read count data from RNA-seq or similar data matrix as input and generates test statistics in consistent format as output. 
</font></div>

# Introduction

Two-group comparison of differential expression (DE) is the most common analysis of transcriptome data. For RNA-seq data, the comparison is usually performed on a gene-level matrix of read counts, with the read counts corresponding to the number of sequencing reads mapped to each gene in each RNA-seq sample. 

Statistical methods that have been applied to two-group DE of RNA-seq data are widely different in terms of their data distribution assumption, input/output format, performance, sensitivity, and user-friendliness. 



<div style="color:darkblue; padding:0 0.1cm">
**Table 1** DE method features. 
</div>

<div align='center', style="padding:0 0.1cm"> 


|Name                                                                                                                          |Call         |Default |Speed  |Paired |Logged |Normalization |Distribution          |Test                                        |Function                 |
|:-----------------------------------------------------------------------------------------------------------------------------|:------------|:-------|:------|:------|:------|:-------------|:---------------------|:-------------------------------------------|:------------------------|
|<a href="https://en.wikipedia.org/wiki/Student%27s_t-test" target="_blank">StudentsT</a>                                      |DeT          |Yes     |Fast   |Yes    |Yes    |No            |Normal                |Student's t test, equal variance            |TDist {stats}            |
|<a href="https://bioconductor.org/packages/release/bioc/html/limma.html" target="_blank">limma</a>                            |DeLimma      |Yes     |Fast   |Yes    |Yes    |No            |Normal                |Empirical Bayes moderation                  |ebayes {limma}           |
|<a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html" target="_blank">edgeR</a>                            |DeEdgeR      |Yes     |Fast   |Yes    |No     |Yes           |Negative binomial     |Exact/Likelihood ratio                      |exactTest {edgeR}        |
|<a href="http://bioconductor.org/packages/release/bioc/html/DESeq.html" target="_blank">DESeq2</a>                            |DeDeSeq      |Yes     |Fast   |Yes    |No     |Yes           |Negative binomial     |Generalized linear model                    |DESeq {DESeq2}           |
|<a href="https://bioconductor.org/packages/release/bioc/html/ABSSeq.html" target="_blank">ABSSeq</a>                          |DeAbsSeq     |No      |Fast   |Yes    |No     |Yes           |Negative binomial     |Sum of counts                               |callDEs {ABSSeq}         |
|<a href="https://www.bioconductor.org/packages/3.3/bioc/html/BGmix.html" target="_blank">BGmix</a>                            |DeBGmix      |No      |Fast   |Yes    |Yes    |No            |Normal                |Bayesian mixture model                      |BGmix {BGmix}            |
|<a href="https://cran.r-project.org/web/packages/PoissonSeq/index.html" target="_blank">PoissonSeq</a>                        |DePoissonSeq |No      |Fast   |Yes    |No     |Yes           |Poisson log-linear    |Poisson goodness-of-fit                     |PS.Main {PoissonSeq}     |
|<a href="http://bioconductor.org/packages/3.3/bioc/html/RBM.html" target="_blank">RBM</a>                                     |DeRBM        |No      |Fast   |No     |Yes    |No            |Normal                |Empirical Bayes & resampling                |RBM_T {RBM}              |
|<a href="http://www.ncbi.nlm.nih.gov/pubmed/24485249" target="_blank">voom</a>                                                |DeVoomLimma  |No      |Fast   |Yes    |No     |Yes           |Log-normal            |Empirical Bayes moderation                  |voom {limma}             |
|<a href="https://en.wikipedia.org/wiki/Welch%27s_t-test" target="_blank">WelchsT</a>                                          |DeWelch      |No      |Fast   |Yes    |Yes    |No            |Normal                |Welch's t test, unequal variance            |TDist {stats}            |
|<a href="http://bioconductor.org/packages/release/bioc/html/DEGseq.html" target="_blank">DEGseq</a>                           |DeDegSeq     |No      |Medium |No     |No     |No            |Binomial/Poisson      |Likelihood Ratio Test                       |DEGexp {DEGseq}          |
|<a href="http://bioconductor.org/packages/release/bioc/html/EBSeq.html" target="_blank">EBSeq</a>                             |DeEbSeq      |No      |Medium |No     |No     |Yes           |Negative Binomial     |Empirical Bayesian                          |EBTest {EBSeq}           |
|<a href="https://bioconductor.org/packages/3.3/bioc/html/NOISeq.html" target="_blank">NOISeq</a>                              |DeNoiSeq     |No      |Medium |No     |No     |Yes           |Nonparametric         |Empirical Bayes                             |noiseqbio {NOISeq}       |
|<a href="https://www.bioconductor.org/packages/release/bioc/html/plgem.html" target="_blank">PLGEM</a>                        |DePlgem      |No      |Medium |Yes    |No     |No            |Normal                |Power Law Global Error Model                |run.plgem {plgem}        |
|<a href="https://bioconductor.org/packages/release/bioc/html/RankProd.html" target="_blank">RankProd</a>                      |DeRankP      |No      |Medium |Yes    |Yes    |No            |Nonparametric         |Rank product                                |RP {RankProd}            |
|<a href="http://statweb.stanford.edu/~tibs/SAM/" target="_blank">SAM</a>                                                      |DeSam        |No      |Medium |Yes    |Yes    |No            |Normal                |Alternative t test with permutation         |samr {samr}              |
|<a href="http://statweb.stanford.edu/~tibs/SAM/" target="_blank">SAMSeq</a>                                                   |DeSamSeq     |No      |Medium |Yes    |No     |No            |Nonparametric         |Wilcoxon with resampling                    |SAMseq {samr}            |
|<a href="http://bioconductor.org/packages/3.3/bioc/html/sSeq.html" target="_blank">sSeq</a>                                   |DeSSeq       |No      |Medium |No     |No     |No            |Negative Binomial     |Shrinkage Approach of Dispersion Estimation |nbTestSH {sSeq}          |
|<a href="https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test" target="_blank">Wilcoxon</a>                                |DeWilcoxon   |No      |Medium |Yes    |Yes    |No            |Nonparametric         |Wilcoxon signed-rank test                   |wilcox.test {stats}      |
|<a href="http://bioconductor.org/packages/release/bioc/html/baySeq.html" target="_blank">baySeq</a>                           |DeBaySeq     |No      |Slow   |Yes    |No     |No            |Negative binomial     |Empirical Bayesian                          |getLikelihoods {baySeq}  |
|<a href="http://bioconductor.org/packages/3.3/bioc/html/bridge.html" target="_blank">bridge</a>                               |DeBridge     |No      |Slow   |No     |Yes    |No            |T/Gaussian            |Bayesian hierarchical model                 |bridge.2samples {bridge} |
|<a href="http://bioconductor.org/packages/release/bioc/html/LMGene.html" target="_blank">LMGene</a>                           |DeLMGene     |No      |Slow   |Yes    |No     |No            |Normal                |Linear model & glog transformation          |genediff {LMGene}        |
|<a href="https://bioconductor.org/packages/release/bioc/html/ALDEx2.html" target="_blank">ALDEx2</a>                          |DeAldex2     |No      |Slower |Yes    |No     |Yes           |Dirichlet             |Welch's t/Wilcoxon/Kruskal Wallace          |aldex {ALDEx2}           |
|<a href="https://www.bioconductor.org/packages/3.3/bioc/html/BADER.html" target="_blank">BADER</a>                            |DeBader      |No      |Slower |No     |No     |Yes           |Overdispersed poisson |Bayesian                                    |BADER {BADER}            |
|<a href="http://bioinformatics.oxfordjournals.org/content/early/2015/04/21/bioinformatics.btv209" target="_blank">edgeRun</a> |DeEdgeRun    |No      |Slower |Yes    |No     |Yes           |Negative binomial     |Exact unconditional                         |UCexactTest {edgeRun}    |
|<a href="http://bioconductor.org/packages/release/bioc/html/tweeDEseq.html" target="_blank">tweeDEseq</a>                     |DeTweeDeSeq  |No      |Slower |No     |No     |Yes           |Poisson-Tweedie       |Poisson-like                                |tweeDE {tweeDEseq}       |


</div>


# Run DE analysis

## Prepare inputs

### Read count matrix

### Grouping samples

### Other parameters

## Run DE analysis online

## Run DE analysis offline

# Browse DE results

## Test statistics 






***
_END OF DOCUMENT_
