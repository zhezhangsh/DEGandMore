\name{DeRNAseq}
\alias{DeRNAseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise DE analysis using RNA-seq read counts
}
\description{
This function uses read count matrix from RNA-seq data and a similar user interface to perform differential expression (DE) analysis between 2 sample groups with one or multiple statistical methods.
}
\usage{
DeRNAseq(ct, grps, paired = FALSE, mthd = 0, min.count = 6, just.stats = TRUE,
         norm.count = c('DESeq', 'TMM', 'RLE', 'QQ', 'UpperQuantile', 'Median', 'TotalCount'),
         norm.logged = c('Loess', 'VST', 'Rlog', 'QQ', 'UpperQuantile', 'Median', force.norm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ct}{
  A matrix of non-negative integers. The columns should be the samples of 2 compared groups and the rows should be a set of genes, exons, or other genomic features. The values of the matrix are numbers of sequencing reads mapped to those features, in terms of RNA-seq data.
  }
  \item{grps}{
  A list of 2 sample name vectors. All sample names should be found within the column names of read count matrix. The vectors correspond to 2 sample groups to be compared. The vectors should be named by group names; otherwise, they will be named 'A' and 'B'.
  }  
  \item{paired}{
  Whether the DE analysis is a paired comparison. If TRUE, the 2 vectors of sample names should have the same length, and in the same order so samples match each other as a pair. Note that some methods don't support paired comparison, so this parameter will be ignored by these methods. 
  } 
  \item{mthd}{
  The method(s) to be used for DE analysis. Call DeMethods() to get a list of available methods. The methods could be named explicitly as a vector of characters, such as c('DeT', 'DeEdgeR'); or it could be a code of a predefined group of methods. Call DeRNAseqMethods() to see full list of available methods.
  
  0: default methods only. Call DeRNAseqMethods(group=0) to see list.
  
  1: fast methods only. Each usually takes less than 10 seconds. Call DeRNAseqMethods(group=1) to see list.
  
  2: fast & medium methods only. Each usually takes less than 1 minute. Call DeRNAseqMethods(group=2) to see list.
  
  3: fast & medium & slow methods. Some methods take up to tens of minutes. Call DeRNAseqMethods(group=3) to see list.
  } 
  \item{min.count}{
  The minimal number of total read counts from all involved samples for a gene to be tested. For genes with less number of reads, their group means, mean difference, and log2(fold change) will be set to 0 while p value and FDR will be set to 1 in the output.
  } 
  \item{num.cluster}{
  The number of clusters to be used for parallel computing. Requires the snow package when it has a value more than 1. Slow methods will be run first. 
  } 
  \item{just.stat}{
  If TRUE, each method will just return the statistical table with the matching rows of the input matrix and 6 columns: group means, mean difference, log2(fold change), p value and FDR; otherwise, all outputs of each method will be returned.
  }
  \item{norm.count}{
  The normalization option for methods use read count data, but don't have their own internal normalization procedure.
  } 
  \item{norm.logged}{
  The normalization option for methods use log2 transformed data as input.
  } 
  \item{force.norm}{
  If TRUE, force normalization of both read counts and logged data no matter whether they are required by the DE methods. 
  } 
}
\details{
  See detailed description of each DE method.
}
\value{
A list of 2 elements:
  \item{inputs}{A list of input data and parameters. The first element is the read count matrix}
  \item{outputs}{A list of outputs from each of the selected DE method. If the _just.stat_ option is TRUE, the output from each method is a 6-column matrix with group means, mean difference, log2(fold change), p value and FDR.}
}
\references{
awsomics.org
}
\author{
Zhe Zhang; zhangz@email.chop.edu
}
\note{
}
\seealso{
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
