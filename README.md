# scAnnot: an automated annotation package for single-cell analysis
![scAnnot logo](https://github.com/eudoraleer/scAnnot/blob/main/man/scAnnot_Logo.png)

__scAnnot__ is an R package created for automated annotation for single-cell data.

## Installation

scAnnot can be installed in R via the command:
```sh
install.packages("devtools")
devtools::install_github("eudoraleer/scAnnot")
```
## Running automated annotation with scAnnot
```sh
library("scAnnot")
result <- scAnnot(project_name = 'scAnnot_Test', ref_degs_path = "ref_degs.csv", data_degs_path = "input_degs")
```
