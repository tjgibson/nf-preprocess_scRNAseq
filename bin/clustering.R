#!/usr/bin/env Rscript

# setup ========================================================================
library(Seurat)
library(anndataR)
library(tidyverse)

# read command line arguments ==================================================
args <- commandArgs(TRUE)
h5ad_fns <- unlist(strsplit(args[1], ","))
out_fn <- args[2]
cluster_npcs <- as.numeric(args[3])
cluster_resolution <- as.numeric(args[4])
integrate_data <- as.logical(args[5])
experiment_name <- args[6]

# read adata and convert to seurat object ======================================
seurat_obj_list <- list()
for (i in seq(h5ad_fns)) {
  adata <- read_h5ad(h5ad_fns[i])
  seurat_obj <- adata$as_Seurat(x_mapping = "counts")
  
  if (i == 1) {
  first_obj <- seurat_obj
} else {
    seurat_obj_list[[i - 1]] <- seurat_obj
}
}

if (length(h5ad_fns > 1)) {
merged_seurat_obj <- merge(first_obj, seurat_obj_list)
} else if (length(h5ad_fns == 1)) {
  merged_seurat_obj <- first_obj
}

rm(first_obj, seurat_obj, seurat_obj_list)

# data normalization, feature selection, scaling and PCA =======================
merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj)
merged_seurat_obj <- ScaleData(merged_seurat_obj)
merged_seurat_obj <- RunPCA(merged_seurat_obj, npcs = cluster_npcs)

# clustering ===================================================================
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, dims = 1:cluster_npcs, reduction = "pca")
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = cluster_resolution, cluster.name = "unintegrated_clusters")

# compute UMAP =================================================================
merged_seurat_obj <- RunUMAP(merged_seurat_obj, dims = 1:cluster_npcs, reduction = "pca", reduction.name = "umap.unintegrated")

# perform RPCA integration =====================================================
if (integrate_data) {
  merged_seurat_obj <- IntegrateLayers(
    object = merged_seurat_obj, method = RPCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.rpca",
    verbose = FALSE
  )
  
  merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "integrated.rpca", dims = 1:cluster_npcs)
  merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = cluster_resolution, cluster.name = "rpca_clusters")
  merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "integrated.rpca", dims = 1:cluster_npcs, reduction.name = "umap.rpca")
}

# write results to file ========================================================
saveRDS(merged_seurat_obj, out_fn)


merged_seurat_obj@meta.data |> 
  as.data.frame() |> 
  rownames_to_column("barcode") |> 
  select(barcode, unintegrated_clusters) |> 
  write_csv(paste0(experiment_name, "res", cluster_resolution,"_unintegrated_clusters", ".csv"))

merged_seurat_obj@reductions$umap.unintegrated@cell.embeddings |> 
  as.data.frame() |> 
  rownames_to_column("barcode") |> 
  write_csv(paste0(experiment_name, "_unintegrated_UMAP.csv"))

if (integrate_data) {
  merged_seurat_obj@meta.data |> 
    as.data.frame() |> 
    rownames_to_column("barcode") |> 
    select(barcode, unintegrated_clusters) |> 
    write_csv(paste0(experiment_name, "res", cluster_resolution,"rpca_clusters", ".csv"))
  
  merged_seurat_obj@reductions$umap.unintegrated@cell.embeddings |> 
    as.data.frame() |> 
    rownames_to_column("barcode") |> 
    write_csv(paste0(experiment_name, "_RPCA_integrated_UMAP.csv"))
}

