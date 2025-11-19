#!/usr/bin/env Rscript

# Setup ========================================================================
suppressPackageStartupMessages({
library(scDblFinder)
library(anndataR)
})

# read command line arguments ==================================================
args <- commandArgs(TRUE)
h5ad_fn <- args[1]
out_fn <- args[2]

# Read in data =================================================================
adata <- read_h5ad(h5ad_fn)

# create sce object ============================================================
sce <- adata$as_SingleCellExperiment(x_mapping = "counts")

# Run scDblFinder ==============================================================
sce <- scDblFinder(sce)

# add doublet score and classification to adata ================================
adata$obs$scDblFinder_score <- sce$scDblFinder.score
adata$obs$scDblFinder_class <- sce$scDblFinder.class
adata$uns$scDblFinder_threshold <- sce@metadata$scDblFinder.threshold

# export adata =================================================================
write_h5ad(adata, out_fn)
