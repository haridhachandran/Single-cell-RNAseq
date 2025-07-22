# ============================
# Single-cell RNA-seq Analysis
# Author: Haridha Sree C
# Date: 2025-07-22
# ============================

# --- Set working directory ---
setwd("/Users/siddhantkalra/Desktop/Decode_Workshop/Sc_RNAseq")

# --- Load libraries ---
library(AnnotationDbi)        # For gene annotation
library(org.Hs.eg.db)         # Human gene annotation database
library(dplyr)                # Data manipulation
library(tidyr)                # Data tidying
library(data.table)           # Data handling
library(Seurat)               # scRNA-seq analysis

# --- Import count data ---
counts <- fread("EXP0001_PCG_beforeQC.txt", sep = "\t", header = TRUE)
counts <- counts[-1, ]        # Remove first row if needed

# Add gene symbols from Ensembl IDs
counts$Gene <- mapIds(org.Hs.eg.db,
                      keys = row.names(counts),
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# Reorder columns to have Gene first
counts <- counts %>% select(Gene, everything()) 

# Remove duplicated genes and rows with missing gene names
counts <- distinct(counts, Gene, .keep_all = TRUE) %>% drop_na(Gene)

# Set gene names as row names and remove Gene column
row.names(counts) <- counts$Gene
counts$Gene <- NULL

# --- Create Seurat object ---
pb <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)

# Calculate percentage of mitochondrial genes
pb[["percent.mt"]] <- PercentageFeatureSet(pb, pattern = "^MT-")

# --- QC Violin plot ---
pdf(file = "Plot1.pdf")
VlnPlot(pb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# --- Normalize data ---
pb <- NormalizeData(pb)

# --- Identify variable features ---
pb <- FindVariableFeatures(pb, selection.method = "vst", nfeatures = 2000)

# --- Visualize top variable genes ---
top10 <- head(VariableFeatures(pb), 10)
plot1 <- VariableFeaturePlot(pb)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# --- Scaling data ---
all.genes <- rownames(pb)
pb <- ScaleData(pb, features = all.genes)

# --- PCA Analysis ---
pb <- RunPCA(pb, features = VariableFeatures(object = pb))
VizDimLoadings(pb, dims = 1:2, reduction = "pca")
DimPlot(pb, reduction = "pca")
DimHeatmap(pb, dims = 1:15, cells = 500, balanced = TRUE)

# --- JackStraw and Elbow plots for dimensionality assessment ---
pb <- JackStraw(pb, num.replicate = 100)
pb <- ScoreJackStraw(pb, dims = 1:20)
JackStrawPlot(pb, dims = 1:15)
ElbowPlot(pb)

# --- Clustering and UMAP visualization ---
pb <- FindNeighbors(pb, dims = 1:10)
pb <- FindClusters(pb, resolution = 0.5)
pb <- RunUMAP(pb, dims = 1:10)
DimPlot(pb, reduction = "umap")

# --- Identify cluster markers ---
pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pb.markers, file = "Markers_info.csv", row.names = FALSE)

# --- Feature visualization for selected genes ---
VlnPlot(pb, features = c("FTL", "RAP1A", "ISG15"))
FeaturePlot(pb, features = c("STATH", "ISG15"), pt.size = 2)

# --- Extract final cluster identities ---
cluster_idents <- Idents(pb)
cluster_df <- data.table::setDT(as.data.frame(cluster_idents), keep.rownames = TRUE)
colnames(cluster_df) <- c("Cell_Id", "Cluster")
