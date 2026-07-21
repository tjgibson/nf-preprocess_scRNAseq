#!/usr/bin/env Rscript

# Setup ========================================================================
suppressPackageStartupMessages({
  library(decontX)
  library(anndataR)
  library(Seurat)
  library(SingleCellExperiment)
})

# read command line arguments ==================================================
args <- commandArgs(TRUE)
filtered_count_fn <- args[1]
raw_count_fn <- args[2]
out_fn <- args[3]
output_corrected_counts <- as.logical(args[4])

# Read in data =================================================================
filtered <- Seurat::Read10X_h5(filtered_count_fn)
raw <- Seurat::Read10X_h5(raw_count_fn)


# create sce object ============================================================
if (is.list(filtered)) {
  sce_mat_list_filtered <- list(counts=filtered$`Gene Expression`)
  sce_mat_list_raw <- list(counts=raw$`Gene Expression`)
} else {
  sce_mat_list_filtered <- list(counts=filtered)
  sce_mat_list_raw <- list(counts=raw)
}

sce <- SingleCellExperiment::SingleCellExperiment(assays = sce_mat_list_filtered)
raw_sce <- SingleCellExperiment::SingleCellExperiment(assays = sce_mat_list_raw)

# Run decontX ==============================================================
sce <- decontX(sce, background = raw_sce)

# creat output adata ==========================================================
if (output_corrected_counts) {
  out_adata <- AnnData(
    X = t(round(sce@assays@data$decontXcounts)), 
    layers = list(uncorrected_counts = t(sce@assays@data$counts))
  )
  
} else {
  out_adata <- AnnData(
    X = t(sce@assays@data$counts), 
    layers = list(decontX_corrected_counts = t(round(sce@assays@data$decontXcounts)))
  )
}

# add doublet score and classification to adata ================================
out_adata$obs$decontX_contamination <- sce$decontX_contamination

# export adata =================================================================
write_h5ad(out_adata, out_fn)
