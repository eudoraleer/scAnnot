% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scAnnot.R
\name{scAnnot}
\alias{scAnnot}
\title{Automated annotation package for single-cell analysis}
\usage{
scAnnot(
  project_name = "scAnnot_Project",
  input_dir = "./",
  output_dir = "./",
  ref_degs_path = NULL,
  ref_degs = NULL,
  ref_degs_format = "general",
  ref_padj_cutoff = 0.05,
  ref_logfc_cutoff = 0.25,
  ref_topn = 30,
  data_degs_path = NULL,
  data_degs = NULL,
  data_degs_format = "general",
  data_padj_cutoff = 0.05,
  result_padj_cutoff = 0.1,
  filter = TRUE
)
}
\arguments{
\item{project_name}{Project name. 'scAnnot_Project' by default.}

\item{input_dir}{Input directory. Default directory
is the current working directory.}

\item{output_dir}{Output directory. Default directory is the
current working directory.}

\item{ref_degs_path}{Path to reference celltype-specific
differentially expressed genes (DEGs). Input format should be CSV or txt.}

\item{ref_degs}{Data frame of the reference celltype-specific degs.}

\item{ref_degs_format}{Data format of the reference celltype-specific degs.
Default is 'general'. Takes input: 'general', 'seurat', 'scanpy'.
'general' format should include data columns: 'group','names','logfoldchanges',
'pvals_adj','pvals', for which 'group' refers to the column for celltype names,
and 'names' refers to the column for gene names.}

\item{ref_padj_cutoff}{P-adjusted threshold. Default is 0.05.}

\item{ref_logfc_cutoff}{Log(fold-change) threshold. Default is 0.25.}

\item{ref_topn}{Top n genes from each reference celltype DEGs to be used
for annotation. Default is 30.}

\item{data_degs_path}{Path to cluster-specific DEGs of the input data to be annotated.
Input format should be CSV or txt.}

\item{data_degs}{Data frame of the cluster-specific degs of the input data.}

\item{data_degs_format}{Data format of the cluster-specific degs of the input data.
Default is 'general'. Takes input: 'general', 'seurat', 'scanpy'.
'general' format should include data columns: 'group','names','logfoldchanges',
'pvals_adj','pvals', for which 'group' refers to the column for celltype names,
and 'names' refers to the column for gene names.}

\item{data_padj_cutoff}{P-adjusted threshold of the cluster-specific degs
of the input data. Default is 0.05.}

\item{result_padj_cutoff}{P-adjusted threshold of the final result.
Default is 0.1.}

\item{filter}{Whether to filter and keep only the top predicted celltype for each cluster.
Default is set to TRUE.}
}
\description{
This function will run automated annotation for single-cell data.
}
