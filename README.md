### Single-cell RNA-seq Data Analysis

This repository contains an R workflow for exploratory analysis and annotation of single-cell RNA-seq data using the Seurat package.

#### Introduction

Single-cell RNA-seq data from human peripheral blood mononuclear cells (PBMCs) of a healthy female donor aged 25-30 were obtained by 10XGenomics from AllCells.

Libraries were generated from ~33,000 cells (23,837 cells recovered) as described in the Chromium Next GEM Single Cell 3' HT Reagent Kits v3.1 User Guide (CG000416 Rev A) using the Chromium X and sequenced on an Illumina NovaSeq 6000 to a read depth of approximately 35,000 mean reads per cell. The data can be downloaded from 10XGenomics database: https://www.10xgenomics.com/datasets/20-k-human-pbm-cs-3-ht-v-3-1-chromium-x-3-1-high-6-1-0

#### Goals

This workflow will perform the following analyses:

-   QC and Filtering 
-   Normalization
-   Identify highly variable features
-   Scale the data
-   Perform Linear dimensionality reduction (PCA)
-   Clustering
-   Non-linear dimensionality reduction (UMAP)
-   Annotation of cell types

#### Documentation and References:

`Seurat`: https://satijalab.org/seurat/

`Seurat` (older versions): https://satijalab.org/seurat/articles/install_v5.html#install-previous-versions-of-seurat

`Seurat` package also requires the `Matrix` as well as `SeuratObject` dependencies. 

They can be installed directly from CRAN and from the developer's website.

install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-1.tar.gz", repos = NULL, type = "source")

remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

Note: `hdf5` command line tool will also need to be installed to process `.HDF5` or `.h5` files.  

First, verify the presence of `hdf5` using the command line console:

h5dump --version

and if necessary, install `hdf5` (e.g. by using homebrew):

brew install hdf5

`SingleR`: https://bioconductor.org/packages/release/bioc/html/SingleR.html

`celldex`: https://bioconductor.org/packages/release/data/experiment/html/celldex.html

`tidyverse`: https://www.tidyverse.org/packages/

`pheatmap`: https://github.com/raivokolde/pheatmap , https://cran.r-project.org/web/packages/pheatmap/index.html


