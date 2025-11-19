# Setup ========================================================================
suppressPackageStartupMessages({
  library(anndataR)
  library(tidyverse)
  library(MASS)
  library(ggnewscale)
  library(plotgardener)
})


# read command line arguments ==================================================
args <- commandArgs(TRUE)
h5ad_fn <- args[1]
plot_fn <- args[2]
min_UMI <-  args[3]
max_UMI <-  args[4]
min_genes_per_cell <- args[5]
max_percent_mito <- args[6]
max_doublet_score <- args[7]
MAD_thresh <- args[8]
percentile_thresh <- args[9]

utils_path <- args[10]

# source utility functions =====================================================
source(utils_path)

# Read in data =================================================================
adata <- read_h5ad(h5ad_fn)

# extract QC data ==============================================================
sc_metadata <- adata$obs |> 
  tibble::rownames_to_column("cell_barcode")

# determine threshold ==========================================================
UMI_thresh_manual <- c(log(min_UMI + 1), log(max_UMI + 1))
UMI_thresh_percentile <- quantile(sc_metadata$log1p_total_counts, c(percentile_thresh, (1 - percentile_thresh)))
UMI_threshold_MAD <- mad_threshold(sc_metadata$log1p_total_counts, MAD_thresh)

mito_thresh_manual <- c(max_percent_mito)
mito_thresh_percentile <- quantile(sc_metadata$pct_counts_mt, (1 - percentile_thresh))
mito_threshold_MAD <- mad_threshold(sc_metadata$pct_counts_mt, MAD_thresh)[2]

gpc_thresh_manual <- log(min_genes_per_cell + 1)
gpc_thresh_percentile <- quantile(sc_metadata$log1p_n_genes_by_counts, percentile_thresh)
gpc_thresh_MAD <- mad_threshold(sc_metadata$log1p_n_genes_by_counts, MAD_thresh)[1]

doublet_thresh_manual <- max_doublet_score
doublet_thresh_class <- adata$uns$scDblFinder_threshold

thresh_data <- tibble(
  metric = c(rep("UMI", 6),rep("n_genes_per_cell", 3), rep("pct_mito", 3), rep("doublet", 2)),
  threshold_strategy = 
    c(
      rep("manual", 2), rep("percentile", 2), rep("MAD", 2), 
      "manual", "percentile", "MAD", 
      "manual", "percentile", "MAD",
      "scDblFinder_threshold", "manual"),
  threshold_value = c(
    UMI_thresh_manual, 
    UMI_thresh_percentile, 
    UMI_threshold_MAD, 
    gpc_thresh_manual,
    gpc_thresh_percentile,
    gpc_thresh_MAD,
    mito_thresh_manual,
    mito_thresh_percentile,
    mito_threshold_MAD,
    doublet_thresh_class,
    doublet_thresh_manual
  )
)

# create blank layout for figure ===============================================
fig_width <-  18
fig_height <- 10

# open pdf
pdf(plot_fn, useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)

# generate plotGardener page
pageCreate(width = fig_width, height = fig_height, default.units = "cm", showGuides = FALSE)

# generate QC plots ============================================================
# UMI violin plot ==============================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5


# get thresholds for plotting
UMI_thresh_data <- thresh_data |> 
  filter(metric == "UMI")

qc_plot <- sc_metadata |> 
  ggplot(aes(x="",y = log1p_total_counts)) + 
  geom_violin(fill = "dodgerblue4") + 
  geom_point(position = position_jitter(seed = 1, width = 0.1), alpha = 0.8, size = 0.05) +
  geom_hline(data = UMI_thresh_data, aes(yintercept = threshold_value, lty = threshold_strategy, color = threshold_strategy)) +
  theme_bw(base_size = 5) +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key.size = unit(2, 'mm'),
    plot.margin=unit(c(0,0,0,0), "mm"), 
    legend.position = "bottom",
  )


# place image on page
plotGG(
  qc_plot,
  x = ref_x + 0.25,
  y = ref_y + 0.25,
  width = 3.5,
  height = 3.5,
  default.units = "cm",
  just = c("left, top")
  
)

# add plot label
label_x <- ref_x + 0.25 + (3.5 / 2)
plotText(
  label = "UMIs per cell", fontsize = 7, fontface = "bold",
  x = label_x, y = ref_y, just = "center", default.units = "cm"
)

# genes per cell violin plot ===================================================
# reference points for positioning figure components
ref_x <- 4.5
ref_y <- 0.5

# get thresholds for plotting
gpc_thresh_data <- thresh_data |> 
  filter(metric == "n_genes_per_cell")

qc_plot <- sc_metadata |> 
  ggplot(aes(x="",y = log1p_n_genes_by_counts)) + 
  geom_violin(fill = "dodgerblue4") + 
  geom_point(position = position_jitter(seed = 1, width = 0.1), alpha = 0.8, size = 0.05) +
  geom_hline(data = gpc_thresh_data, aes(yintercept = threshold_value, lty = threshold_strategy, color = threshold_strategy)) +
  theme_bw(base_size = 5) +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key.size = unit(2, 'mm'),
    plot.margin=unit(c(0,0,0,0), "mm"), 
    legend.position = "bottom",
    legend.title = element_blank()
  )


# place image on page
plotGG(
  qc_plot,
  x = ref_x + 0.25,
  y = ref_y + 0.25,
  width = 3.5,
  height = 3.5,
  default.units = "cm",
  just = c("left, top")
  
)

# add plot label
label_x <- ref_x + 0.25 + (3.5 / 2)
plotText(
  label = "Genes per cell", fontsize = 7, fontface = "bold",
  x = label_x, y = ref_y, just = "center", default.units = "cm"
)

# percent mitochondira violin plot =============================================
# reference points for positioning figure components
ref_x <- 8.5
ref_y <- 0.5

mito_thresh_data <- thresh_data |> 
  filter(metric == "pct_mito")


qc_plot <- sc_metadata |> 
  ggplot(aes(x="",y = pct_counts_mt)) + 
  geom_violin(fill = "dodgerblue4") + 
  geom_point(position = position_jitter(seed = 1, width = 0.1), alpha = 0.8, size = 0.05) +
  geom_hline(data = mito_thresh_data, aes(yintercept = threshold_value, lty = threshold_strategy, color = threshold_strategy)) +
  theme_bw(base_size = 5) +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key.size = unit(2, 'mm'),
    plot.margin=unit(c(0,0,0,0), "mm"), 
    legend.position = "bottom",
    legend.title = element_blank()
  )

# place image on page
plotGG(
  qc_plot,
  x = ref_x + 0.25,
  y = ref_y + 0.25,
  width = 3.5,
  height = 3.5,
  default.units = "cm",
  just = c("left, top")
  
)

# add plot label
label_x <- ref_x + 0.25 + (3.5 / 2)
plotText(
  label = "Percent mitochondrial", fontsize = 7, fontface = "bold",
  x = label_x, y = ref_y, just = "center", default.units = "cm"
)

# doublet score violin plot ====================================================
# reference points for positioning figure components
ref_x <- 12.5
ref_y <- 0.5

# get thresholds for plotting
doublet_thresh_data <- thresh_data |> 
  filter(metric == "doublet")

qc_plot <- sc_metadata |> 
  ggplot(aes(x="",y = scDblFinder_score)) + 
  geom_violin(fill = "dodgerblue4") + 
  geom_point(position = position_jitter(seed = 1, width = 0.1), alpha = 0.8, size = 0.05) +
  geom_hline(data = doublet_thresh_data, aes(yintercept = threshold_value, lty = threshold_strategy, color = threshold_strategy)) +
  theme_bw(base_size = 5) +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key.size = unit(2, 'mm'),
    plot.margin=unit(c(0,0,0,0), "mm"), 
    legend.position = "bottom",
    legend.title = element_blank()
  )

# place image on page
plotGG(
  qc_plot,
  x = ref_x + 0.25,
  y = ref_y + 0.25,
  width = 3.5,
  height = 3.5,
  default.units = "cm",
  just = c("left, top")
  
)

# add plot label
label_x <- ref_x + 0.25 + (3.5 / 2)
plotText(
  label = "Doublet score (scDblFinder)", fontsize = 7, fontface = "bold",
  x = label_x, y = ref_y, just = "center", default.units = "cm"
)


# counts vs genes scatter w MAD thresholds =====================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 4.5

# compute point density for plotting
p_data <- sc_metadata
p_data$density <- get_density(p_data$log1p_total_counts, p_data$log1p_n_genes_by_counts, n = 500)

# get MAD threshold data
UMI_MAD_thresholds <- tibble(
  metric = rep("UMI", 4),
  threshold_strategy = 
    c(
      rep("3X MAD", 2), rep("5X MAD", 2)),
  threshold_value = c(
    mad_threshold(sc_metadata$log1p_total_counts, 3), 
    mad_threshold(sc_metadata$log1p_total_counts, 5) 
  )
)

gpc_MAD_thresholds <- tibble(
  metric = rep("genes_per_cell", 2),
  threshold_strategy = c("3X MAD","5X MAD"),
  threshold_value = c(
    mad_threshold(sc_metadata$log1p_n_genes_by_counts, 3)[1], 
    mad_threshold(sc_metadata$log1p_n_genes_by_counts, 5)[1] 
  )
)

# generate plot
qc_plot <- p_data |> 
  ggplot(aes(x = log1p_total_counts ,y = log1p_n_genes_by_counts, color = density)) + 
  geom_point(size = 0.01) +
 scale_color_viridis_c(guide = "none") +
  new_scale_color() +
  geom_vline(data = UMI_MAD_thresholds, aes(xintercept = threshold_value, lty = threshold_strategy, color = threshold_strategy)) +
  geom_hline(data = gpc_MAD_thresholds, aes(yintercept = threshold_value, lty = threshold_strategy, color = threshold_strategy)) +
  theme_bw(base_size = 5) +
  theme(
    legend.key.size = unit(2, 'mm'),
    plot.margin=unit(c(0,0,0,0), "mm"), 
    legend.position = "bottom",
    legend.title = element_blank()
  )

# place image on page
plotGG(
  qc_plot,
  x = ref_x + 0.25,
  y = ref_y + 0.25,
  width = 3.5,
  height = 3.5,
  default.units = "cm",
  just = c("left, top")
  
)

# add plot label
label_x <- ref_x + 0.25 + (3.5 / 2)
plotText(
  label = "MAD thresholds", fontsize = 7, fontface = "bold",
  x = label_x, y = ref_y, just = "center", default.units = "cm"
)



# counts vs genes scatter w percetile thresholds ===============================
# reference points for positioning figure components
ref_x <- 4.5
ref_y <- 4.5

# get percentile threshold data
UMI_percentile_thresholds <- tibble(
  metric = rep("UMI", 6),
  threshold_strategy = 
    c(
      rep("95th percentile", 2), rep("99th percentile", 2), rep("99.9th percentile", 2)),
  threshold_value = c(
    quantile(sc_metadata$log1p_total_counts, c(0.05, 0.95)), 
    quantile(sc_metadata$log1p_total_counts, c(0.01, 0.99)),
    quantile(sc_metadata$log1p_total_counts, c(0.001, 0.999)) 
  )
)

gpc_percentile_thresholds <- tibble(
  metric = rep("genes_per_cell", 3),
  threshold_strategy = 
    c("95th percentile", "99th percentile", "99.9th percentile"),
  threshold_value = c(
    quantile(sc_metadata$log1p_n_genes_by_counts, c(0.05)), 
    quantile(sc_metadata$log1p_n_genes_by_counts, c(0.01)),
    quantile(sc_metadata$log1p_n_genes_by_counts, c(0.001)) 
  )
)

# generate plot
qc_plot <- p_data |> 
  ggplot(aes(x = log1p_total_counts ,y = log1p_n_genes_by_counts, color = density)) + 
  geom_point(size = 0.01) +
  scale_color_viridis_c(guide = "none") +
  new_scale_color() +
  geom_vline(data = UMI_percentile_thresholds, aes(xintercept = threshold_value, lty = threshold_strategy, color = threshold_strategy)) +
  geom_hline(data = gpc_percentile_thresholds, aes(yintercept = threshold_value, lty = threshold_strategy, color = threshold_strategy)) +
  theme_bw(base_size = 5) +
  theme(
    legend.key.size = unit(2, 'mm'),
    plot.margin=unit(c(0,0,0,0), "mm"), 
    legend.position = "bottom",
    legend.title = element_blank()
  )

# place image on page
plotGG(
  qc_plot,
  x = ref_x + 0.25,
  y = ref_y + 0.25,
  width = 3.5,
  height = 3.5,
  default.units = "cm",
  just = c("left, top")
  
)

# add plot label
label_x <- ref_x + 0.25 + (3.5 / 2)
plotText(
  label = "Percentile thresholds", fontsize = 7, fontface = "bold",
  x = label_x, y = ref_y, just = "center", default.units = "cm"
)

# counts vs genes scatter w percent mitochonrial ===============================
# reference points for positioning figure components
ref_x <- 8.5
ref_y <- 4.5



qc_plot <- sc_metadata |> 
  ggplot(aes(x = log1p_total_counts ,y = log1p_n_genes_by_counts, color = pct_counts_mt)) + 
  geom_point(size = 0.01) +
  scale_color_viridis_c() +
  theme_bw(base_size = 5) +
  theme(
    legend.key.size = unit(2, 'mm'),
    plot.margin=unit(c(0,0,0,0), "mm"), 
    legend.position = "bottom",
    legend.title = element_blank()
  )

# place image on page
plotGG(
  qc_plot,
  x = ref_x + 0.25,
  y = ref_y + 0.25,
  width = 3.5,
  height = 3.5,
  default.units = "cm",
  just = c("left, top")
  
)

# add plot label
label_x <- ref_x + 0.25 + (3.5 / 2)
plotText(
  label = "Percent Mitochondrial", fontsize = 7, fontface = "bold",
  x = label_x, y = ref_y, just = "center", default.units = "cm"
)


# counts vs genes scatter w doublet score  =====================================
# reference points for positioning figure components
ref_x <- 12.5
ref_y <- 4.5



qc_plot <- 
  sc_metadata |> 
  ggplot(aes(x = log1p_total_counts ,y = log1p_n_genes_by_counts, color = scDblFinder_score)) + 
  geom_point(size = 0.05) +
  scale_color_viridis_c()  +
  theme_bw(base_size = 5) +
  theme(
    legend.key.size = unit(2, 'mm'),
    plot.margin=unit(c(0,0,0,0), "mm"), 
    legend.position = "bottom",
    legend.title = element_blank()
  )

# place image on page
plotGG(
  qc_plot,
  x = ref_x + 0.25,
  y = ref_y + 0.25,
  width = 3.5,
  height = 3.5,
  default.units = "cm",
  just = c("left, top")
  
)

# add plot label
label_x <- ref_x + 0.25 + (3.5 / 2)
plotText(
  label = "Doublet score (scDblFinder)", fontsize = 7, fontface = "bold",
  x = label_x, y = ref_y, just = "center", default.units = "cm"
)

# close graphics device ========================================================
dev.off()

