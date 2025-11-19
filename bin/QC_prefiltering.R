#!/usr/bin/env Rscript

# Setup ========================================================================
suppressPackageStartupMessages({
  library(anndataR)
  library(tidyverse)
})

# read command line arguments ==================================================
args <- commandArgs(TRUE)
h5ad_fn <- args[1]
out_fn <- args[2]
UMI_filter_strategy<- args[3]
gene_per_cell_filter_strategy <- args[4]
mito_filter_strategy <- args[5]
doublet_filter_strategy <- args[6]
min_UMI <- args[7]
max_UMI <- args[8]
min_genes_per_cell <- args[9]
max_percent_mito <- args[10]
max_doublet_score <- args[11]
MAD_thresh <- args[12]
percentile_thresh <- args[13]


utils_path <- args[14]

filter_stats_fn <- args[15]

# source utility functions =====================================================
source(utils_path)


# Read in data =================================================================
adata <- read_h5ad(h5ad_fn)

# define filter thresholds =====================================================
# intialize list to store thresholds
thresholds <- list(
  UMI_threshold = list(lower = NA, upper = NA),
  gpc_threshold = list(lower = NA),
  mito_threshold = list(upper = NA),
  doublet_threshold = list(upper = NA)
)

# check for valid filter strategy params
if (!(UMI_filter_strategy %in% c("MAD", "percentile", "manual"))) {
  stop("UMI_filter_strategy should be one of 'MAD', 'percentile', or 'manual'")
}

if (!(gene_per_cell_filter_strategy %in% c("MAD", "percentile", "manual"))) {
  stop("gene_per_cell_filter_strategy should be one of 'MAD', 'percentile', or 'manual'")
}

if (!(mito_filter_strategy %in% c("MAD", "percentile", "manual"))) {
  stop("mito_filter_strategy should be one of 'MAD', 'percentile', or 'manual'")
}

if (!(doublet_filter_strategy %in% c("automatic","manual"))) {
  stop("mito_filter_strategy should be one of 'automatic' or 'manual'")
}

# set UMI threshold values based on input params
if (UMI_filter_strategy == "MAD") {
  thresholds$UMI_threshold$lower <- mad_threshold(adata$obs$log1p_total_counts, MAD_thresh)[1]
  thresholds$UMI_threshold$upper <- mad_threshold(adata$obs$log1p_total_counts, MAD_thresh)[2]
} else if (UMI_filter_strategy == "percentile") {
  thresholds$UMI_threshold$lower <- quantile(adata$obs$log1p_total_counts, percentile_thresh)
  thresholds$UMI_threshold$upper <- quantile(adata$obs$log1p_total_counts, (1 - percentile_thresh))
} else if (UMI_filter_strategy == "manual") {
  thresholds$UMI_threshold$lower <- log(min_UMI + 1)
  thresholds$UMI_threshold$upper <- log(max_UMI + 1)
}

# set genes per cell threshold values based on input params
if (gene_per_cell_filter_strategy == "MAD") {
  thresholds$gpc_threshold$lower <- mad_threshold(adata$obs$log1p_n_genes_by_counts, MAD_thresh)[1]
} else if (gene_per_cell_filter_strategy == "percentile") {
  thresholds$gpc_threshold$lower <- quantile(adata$obs$log1p_n_genes_by_counts, percentile_thresh)
} else if (gene_per_cell_filter_strategy == "manual") {
  thresholds$gpc_threshold$lower <- log(min_genes_per_cell + 1)
}

# set percent mitochondrial threshold values based on input params
if (mito_filter_strategy == "MAD") {
  thresholds$mito_threshold$upper <- mad_threshold(adata$obs$pct_counts_mt, MAD_thresh)[2]
} else if (mito_filter_strategy == "percentile") {
  thresholds$mito_threshold$upper <- quantile(adata$obs$pct_counts_mt, (1 - percentile_thresh))
} else if (mito_filter_strategy == "manual") {
  thresholds$mito_threshold$upper <- max_percent_mito
}

# set doublet score threshold values based on input params
if (doublet_filter_strategy == "automatic") {
  thresholds$doublet_threshold$upper <- adata$uns$scDblFinder_threshold
} else if (mito_filter_strategy == "manual") {
  thresholds$doublet_threshold$upper <- max_doublet_score
}

# perform filtering ============================================================
# initialize table to store filtering stats
filter_stats <- tibble(
  filter_stage = c("total cells before filtering", "low_UMIs", "high_UMIs", "low_genes_per_cell", "high_percent_mitochondrial", "high_doublet_score", "total cells after filtering"),

    filter_threshold = c(
      NA, 
      exp(thresholds$UMI_threshold$lower) - 1, 
      exp(thresholds$UMI_threshold$upper) - 1, 
      exp(thresholds$gpc_threshold$lower) - 1, 
      thresholds$mito_threshold$upper, 
      thresholds$doublet_threshold$upper,
      NA
      ),

    n_cells_removed = as.integer(rep(NA, times = 7)),
  n_cells_remaining = as.integer(rep(NA, times = 7))
  
)

filter_stats$n_cells_remaining[1] <- adata$n_obs()

# perform UMI filtering
keep_cells <- adata$obs$log1p_total_counts > thresholds$UMI_threshold$lower
filtered_adata <- adata[keep_cells]

filter_stats$n_cells_remaining[2] <- filtered_adata$n_obs()
filter_stats$n_cells_removed[2] <- sum(!keep_cells)

keep_cells <- filtered_adata$obs$log1p_total_counts < thresholds$UMI_threshold$upper
filtered_adata <- filtered_adata[keep_cells]

filter_stats$n_cells_remaining[3] <- filtered_adata$n_obs()
filter_stats$n_cells_removed[3] <- sum(!keep_cells)

# perform genes per cell filtering
keep_cells <- filtered_adata$obs$log1p_n_genes_by_counts > thresholds$gpc_threshold$lower
filtered_adata <- filtered_adata[keep_cells]

filter_stats$n_cells_remaining[4] <- filtered_adata$n_obs()
filter_stats$n_cells_removed[4] <- sum(!keep_cells)

# perform mitochondrial filtering
keep_cells <- filtered_adata$obs$pct_counts_mt < thresholds$mito_threshold$upper
filtered_adata <- filtered_adata[keep_cells]

filter_stats$n_cells_remaining[5] <- filtered_adata$n_obs()
filter_stats$n_cells_removed[5] <- sum(!keep_cells)

# perform doublet filtering
keep_cells <- filtered_adata$obs$scDblFinder_score < thresholds$doublet_threshold$upper
filtered_adata <- filtered_adata[keep_cells]

filter_stats$n_cells_remaining[6] <- filtered_adata$n_obs()
filter_stats$n_cells_removed[6] <- sum(!keep_cells)

# add total cells after filtering
filter_stats$n_cells_remaining[7] <- filtered_adata$n_obs()
filter_stats$n_cells_removed[7] <- adata$n_obs() - filtered_adata$n_obs()


# write filtering stats to file ================================================
write_csv(filter_stats, filter_stats_fn)

# write filtered adata to file =================================================
write_h5ad(filtered_adata, out_fn)

