#!/usr/bin/env python3

# import libraries =========================================================================
import scanpy as sc
import anndata as ad
import sys

# parse command line arguments =============================================================
cellranger_h5_path = sys.argv[1]
out_fn = sys.argv[2]
mito_prefix = sys.argv[3]

## create spatialdata object ===============================================================
print("importing cellranger count table")
adata = sc.read_10x_h5(cellranger_h5_path)
    
# make the var names unique ================================================================
print("making variable names unique")
adata.var_names_make_unique()

# compute basic QC metrics =================================================================
    # mitochondrial genes, "MT-" for human, "mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)    
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)

# save h5ad file ===========================================================================
adata.write_h5ad(out_fn)