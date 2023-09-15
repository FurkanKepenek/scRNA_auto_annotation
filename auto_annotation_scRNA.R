# Script for automatic cell annotation in single cell RNA-seq
# https://www.10xgenomics.com/resources/datasets/20-k-human-pbm-cs-3-ht-v-3-1-chromium-x-3-1-high-6-1-0 /// Data is available here. 
# There could be a form to reach data, fill the form then download the "Feature / cell matrix HDF5 (filtered)" data
# needed packages will be installed from bioconductor repo with codes given below

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("SingleR")


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("celldex")

library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(viridis)
library(pheatmap)

# Human peripheral blood mononuclear cells (PBMCs) of a healthy female donor aged 25-30 were obtained by 10x Genomics from AllCells.

# read hdf5 object

hdf5_obj <- Read10X_h5(filename = "20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
                       use.names = TRUE,
                       unique.features = TRUE)

# create seurat object of hdf5 object

pmbc_seurat <- CreateSeuratObject(counts = hdf5_obj)

# calculate the mitochondrial percentage of data

pmbc_seurat$mitoPercent <- PercentageFeatureSet(pmbc_seurat, pattern = '^MT-')

# filter data according to given parameters below
pmbc_seurat_filtered <- subset(pmbc_seurat, subset = nCount_RNA > 800 &
         nFeature_RNA > 500 &
         mitoPercent < 10)


# --------------------pre processing of seurat object---------------------------

# normalize the data 

pmbc_seurat_filtered <- NormalizeData(object = pmbc_seurat_filtered)

# find variable features 

pmbc_seurat_filtered <- FindVariableFeatures(object = pmbc_seurat_filtered)

# scale the data

pmbc_seurat_filtered <- ScaleData(object = pmbc_seurat_filtered)

# start PCA analysis

pmbc_seurat_filtered <- RunPCA(object = pmbc_seurat_filtered)

# find neighbors according to PCA analysis

pmbc_seurat_filtered <- FindNeighbors(object = pmbc_seurat_filtered, dims = 1:20)

# find clusters according to neighbors 

pmbc_seurat_filtered <- FindClusters(object = pmbc_seurat_filtered)

# start UMAP analysis

pmbc_seurat_filtered <- RunUMAP(object = pmbc_seurat_filtered, dims = 1:20)

# -----------------------------Visualization------------------------------------

# have a look at final seurat data

View(pmbc_seurat_filtered@meta.data)

# create plot to visualize seurat data, reduction is "umap" because before UMAP analysis used
DimPlot(pmbc_seurat_filtered, reduction = 'umap')

# --------------------------------Annotation Part-------------------------------

# create reference variable for annotation, there are different datasets for 
# needs for this data Human Primary Cell Atlas Data will be usefull
# thanks to "celldex" package we are reaching these datasets with one line of code

ref <- celldex::HumanPrimaryCellAtlasData()

# observe reference data

View(as.data.frame(colData(ref)))



# get count slot from filtered data

pbmc_counts <- GetAssayData(pmbc_seurat_filtered, slot = 'counts')


# run SingleR function for automatic cell annotation, test will be our count data
# ref will be our reference which taken from "celldex", labels will be the lables 
# from reference data -there will be label.main column- 

single_data <- SingleR(test = pbmc_counts,
        ref = ref,
        labels = ref$label.main)

# observe the result, SingleR gives S4 object as a result

single_data

# provide labels for main data

pmbc_seurat_filtered$singleR.labels <- single_data$labels[match(rownames(pmbc_seurat_filtered@meta.data), rownames(single_data))]

# create plot with labelled cells  

DimPlot(pmbc_seurat_filtered, reduction = 'umap', group.by = 'singleR.labels')

# see the scores of cells

single_data
single_data$scores

# to see ambiguity create score heatmap, each little part represents cells 

plotScoreHeatmap(single_data)


# according to delta scores between cells create delta distribution plot

plotDeltaDistribution(single_data)

# comparing cell scores with unsupervised clustering

tab <- table(Assigned=single_data$labels, Clusters=pmbc_seurat_filtered$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))




