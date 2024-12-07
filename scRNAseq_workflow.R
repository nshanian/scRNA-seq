### Exploratory Analysis and Annotation of single-cell RNA-seq data

### Introduction

This R workflow will perform exploratory analysis and annotation of single-cell RNA-seq data using the `Seurat`, `celldex` and `SingleR` packages. 

Human peripheral blood mononuclear cells (PBMCs) of a healthy female donor aged 25-30 were obtained by 10XGenomics from AllCells.

Libraries were generated from ~33,000 cells (23,837 cells recovered) as described in the Chromium Next GEM Single Cell 3' HT Reagent Kits v3.1 User Guide (CG000416 Rev A) using the Chromium X and sequenced on an Illumina NovaSeq 6000 to a read depth of approximately 35,000 mean reads per cell. The data can be downloaded from 10XGenomics database:

https://www.10xgenomics.com/datasets/20-k-human-pbm-cs-3-ht-v-3-1-chromium-x-3-1-high-6-1-0

### Goals

This workflow will perform the following analyses:

-   QC and Filtering 
-   Normalization
-   Identify highly variable features
-   Scale the data
-   Perform Linear dimensionality reduction (PCA)
-   Clustering
-   Non-linear dimensionality reduction (UMAP)
-   Annotation of cell types

### Setup

As in a regular R script in RStudio, a single line of code can be run with Command-Enter (Mac OS) or Ctrl-Enter (Windows). Whole chunks of code can be run with Command(/Ctrl) + Shift + Enter **or** by clicking the green "\>" button in the top-right corner of the chunk. Alternatively, these options can also be implemented by selecting lines of code and choosing the desired option from the drop-down menu in the "Run" tab, in the top-right corner of the of Source section of RStudio.

To comment or uncomment a series of lines, highlight the lines and use Command(/Ctrl) + Shift + C.

```{css, echo = FALSE}
pre, code {white-space:pre !important; overflow-x:auto}
```


```{r setup, include=FALSE}
# when you "knit" this file, do you want the resulting PDF to print the code in each chunk (TRUE = yes)?
knitr::opts_chunk$set(echo = TRUE)

################################################################################
# set your working directory
####CHANGE THIS TO THE APPROPRIATE PATH
knitr::opts_knit$set(root.dir = '~/Desktop/10XGenomics/')
################################################################################
# note that outside of an Rmd code chunk, use `setwd()` to set the working directory in R
```

R packages `SingleR`, `celldex`, `Seurat`, `tidyverse` and `pheatmap` are required for this workflow. 

They can be installed using the following commands:

```{r install packages, eval=F}
#install.packages("SingleR")
#install.packages("celldex")
#install.packages("Seurat")
#install.packages("tidyverse")
#install.packages("pheatmap")

# If install.packages("packagename") command fails in newer versions of R uncomment and run the commands:
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("packagename")
```

`Seurat` package also requires the `Matrix` as well as `SeuratObject` dependencies. 

They can be installed directly from CRAN and from the developer's website.

```{r install dependencies}
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-1.tar.gz", repos = NULL, type = "source")

# remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

# Note: `hdf5` command line tool will also need to be installed to process `.HDF5` or `.h5` files.  

# First, verify the presence of `hdf5` using the command line console and if necessary, install `hdf5`, e.g. by using homebrew:

# h5dump --version

# brew install hdf5
```

Once installed, load the packages.

```{r load libraries}
library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)
```


```{r import data}
# Import data from 10X CellRanger in .HDF5 format 
hdf5_obj <- Read10X_h5(filename = '~/Desktop/10XGenomics/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5',
                       use.names = TRUE,
                       unique.features = TRUE)
pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)
```

### QC and Filtering

##### It's good practice to filter out cells where genes identified and genes with expression across cells fall below a certain threshold.

```{r filter data}
# QC and Filtering
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
         nFeature_RNA > 500 &
         mitoPercent < 10)
```

### Normalization

```{r normalize data}
#pbmc.seurat.filtered <- NormalizeData(pbmc.seurat.filtered, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
pbmc.seurat.filtered <- NormalizeData(pbmc.seurat.filtered)
str(pbmc.seurat.filtered)
```

### Identify highly variable features

```{r identify variable features}
# Identify highly variable features
pbmc.seurat.filtered <- FindVariableFeatures(pbmc.seurat.filtered, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc.seurat.filtered), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc.seurat.filtered)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
```
### Scale the data

```{r scaling the data}
# Scaling the data
all.genes <- rownames(pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(pbmc.seurat.filtered, features = all.genes)

str(pbmc.seurat.filtered)
```
### Perform Linear dimensionality reduction (PCA)

```{r linear dimensionality reduction}
# Perform Linear dimensionality reduction 
pbmc.seurat.filtered <- RunPCA(pbmc.seurat.filtered, features = VariableFeatures(object = pbmc.seurat.filtered))

# visualize PCA results
print(pbmc.seurat.filtered[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(pbmc.seurat.filtered, dims = 1, cells = 500, balanced = TRUE)
```

#### Creating an elbow or Scree plot is very useful for selecting the optimum number of dimensions or principal components.

```{r determine dimensionality of the data}
# determine dimensionality of the data
ElbowPlot(pbmc.seurat.filtered)
```

### Clustering

```{r clustering}
# Clustering
# In this case standard deviation or variance no longer explains the components or dimensions greater than 15
pbmc.seurat.filtered <- FindNeighbors(pbmc.seurat.filtered, dims = 1:15)

pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
```

### Non-linear dimensionality reduction (UMAP)

```{r non-linear dimensionality reduction}
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:15)
```

```{r preprocessing}
# Alternatively, all preprocessing steps can be run as one chunk

#pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
#pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
#pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
#pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
#pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:15)
#pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
#pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:15)
```

```{r visualize the results}
#View(pbmc.seurat.filtered@meta.data)
DimPlot(pbmc.seurat.filtered, reduction = 'umap')
```

### Annotation 

#### Unlike marker-based annotation that relies on traditional information of canonical marker genes that are specifically expressed in known cell types, this workflow uses the reference-based annotation strategy, which utilizes comprehensive gene expression profiles of expertly annotated reference datasets, like the Human Cell Atlas.

#### Obtain reference data using `celldex` package

```{r annotation reference}
# obtain reference data
# expression values are provided as log normalized counts
ref <- celldex::HumanPrimaryCellAtlasData()
#View(as.data.frame(colData(ref)))
```

#### Annotate cells using `SingleR` package

```{r annotation}
# run SingleR (default mode)
# by default SingleR will perform annotation of each individual cell in the test dataset

pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')

pred <- SingleR(test = pbmc_counts,
        ref = ref,
        labels = ref$label.main)
pred
```

#### Label annotation results

```{r annotation lables}
pbmc.seurat.filtered$singleR.labels <- pred$labels[match(rownames(pbmc.seurat.filtered@meta.data), rownames(pred))]
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'singleR.labels')
```

#### Run annotation dignostics

```{r annotation diganostics}
# Annotation diagnostics
# Based on the scores within cells
pred
#pred$scores

plotScoreHeatmap(pred)

# Based on deltas across cells 
plotDeltaDistribution(pred)
```

#### Compare annotation results

```{r compare annotation results}
# Compare annotation results to unsupervised clustering

tab <- table(Assigned=pred$labels, Clusters=pbmc.seurat.filtered$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))
```

