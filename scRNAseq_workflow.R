# This R script will perform exploratory analysis and annotation of single-cell RNA-seq data
# Source: Human peripheral blood mononuclear cells (PBMCs) of a healthy female donor aged 25-30 were obtained by 10XGenomics from AllCells.

setwd("~/Desktop/10XGenomics/")

# Load the necessary libraries

library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

# Import data from 10X CellRanger in .HDF5 format 

hdf5_obj <- Read10X_h5(filename = '~/Desktop/10XGenomics/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5',
                       use.names = TRUE,
                       unique.features = TRUE)
pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)

# QC and Filtering

pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 mitoPercent < 10)

# Normalize the data

#pbmc.seurat.filtered <- NormalizeData(pbmc.seurat.filtered, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
pbmc.seurat.filtered <- NormalizeData(pbmc.seurat.filtered)
str(pbmc.seurat.filtered)

# Identify highly variable features
pbmc.seurat.filtered <- FindVariableFeatures(pbmc.seurat.filtered, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc.seurat.filtered), 10)


# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc.seurat.filtered)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


# Scaling the data
all.genes <- rownames(pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(pbmc.seurat.filtered, features = all.genes)

str(pbmc.seurat.filtered)

# Creating an elbow or Scree plot is very useful for selecting the optimum number of dimensions or principal components.

# Determine dimensionality of the data
ElbowPlot(pbmc.seurat.filtered)

# Clustering

# In this case standard deviation or variance no longer explains the components or dimensions greater than 15
pbmc.seurat.filtered <- FindNeighbors(pbmc.seurat.filtered, dims = 1:15)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)


# Non-linear dimensionality reduction}
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:15)


# Visualize the results
#View(pbmc.seurat.filtered@meta.data)
DimPlot(pbmc.seurat.filtered, reduction = 'umap')


# Reference-based annotation
# obtain reference data
# expression values are provided as log normalized counts
ref <- celldex::HumanPrimaryCellAtlasData()
#View(as.data.frame(colData(ref)))


# run SingleR (default mode)
# by default SingleR will perform annotation of each individual cell in the test dataset

pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')

pred <- SingleR(test = pbmc_counts,
                ref = ref,
                labels = ref$label.main)
pred

# Label Annotation results
pbmc.seurat.filtered$singleR.labels <- pred$labels[match(rownames(pbmc.seurat.filtered@meta.data), rownames(pred))]
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'singleR.labels')


# Run annotation diagnostics
# Based on the scores within cells
pred
#pred$scores

plotScoreHeatmap(pred)

# Based on deltas across cells 
plotDeltaDistribution(pred)


# Compare annotation results to unsupervised clustering
tab <- table(Assigned=pred$labels, Clusters=pbmc.seurat.filtered$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))


