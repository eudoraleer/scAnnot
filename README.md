# scAnnot: an automated annotation package for single-cell analysis
![scAnnot logo](https://github.com/eudoraleer/scAnnot/blob/main/man/scAnnot_Logo.png)

__scAnnot__ is an R package created for automated annotation for single-cell data.

## Installation

scAnnot can be installed in R via the command:
```sh
install.packages("devtools")
devtools::install_github("eudoraleer/scAnnot")
```
## Preparation of input files
scAnnot takes in the following input files:
1. Differentially expressed genes (DEGs) of a set of reference cell types (granularity depends on the number of clusters to be annotated). The file could be a direct output file from Scanpy or Seurat. Input parameters: ref_degs/ref_degs_path (see ?scAnnot for more help.)
2. Differentially expressed genes (DEGs) of the clusters to be annotated. The file could be a direct output file from Scanpy or Seurat. Input parameters: data_degs/data_degs_path (see ?scAnnot for more help.)

## Running automated annotation with scAnnot
```sh
library("scAnnot")
?scAnnot
result <- scAnnot(project_name = 'scAnnot_Test', ref_degs_path = "ref_degs.csv", data_degs_path = "input_degs.csv")
```
