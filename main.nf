#! /usr/bin/env nextflow

log.info """
	scRNAseq preprocessing pipeline
	===================================
	samplesheet: ${params.samplesheet}
	results_dir: ${params.results_dir}

	UMI_filter_strategy: ${params.UMI_filter_strategy}
	gene_per_cell_filter_strategy: ${params.gene_per_cell_filter_strategy}
	mito_filter_strategy: ${params.mito_filter_strategy}
	doublet_filter_strategy: ${params.doublet_filter_strategy}
	min_UMI: ${params.min_UMI}
	max_UMI: ${params.max_UMI}
	min_genes_per_cell: ${params.min_genes_per_cell}
	max_percent_mito: ${params.max_percent_mito} 
	max_doublet_score: ${params.max_doublet_score} 
	MAD_thresh: ${params.MAD_thresh} 
	percentile_thresh: ${params.percentile_thresh} 
	batch_correction_method: ${params.batch_correction_method} 
	mitochondrial_prefix: ${params.mitochondrial_prefix} 
	R_utils_file: ${params.R_utils_file} 
	cluster_npcs: ${params.cluster_npcs} 
  	cluster_resolution: ${params.cluster_resolution} 
  	integrate_samples: ${params.integrate_samples} 
	ambient RNA correction: ${params.correct_ambient_RNA} 
	"""
	.stripIndent()

/* input files:	
 * samplesheet
 * cellranger h5 file
 */
 


/*
 * estimate ambient RNA contamination
 */

 process decontX {
	tag "$meta.sample"
	container = "tjmgibson/scrnaseq_preprocess:v2"
	publishDir "${params.results_dir}/h5ad_files/unfiltered/", mode: 'copy'
	
	input:
	tuple val(meta), path(filtered_counts), path(raw_counts)
	
	output:
	tuple val(meta), path("${meta.sample}_decontx.h5ad")
	
	script:
    """
	decontX.R ${h5ad} ${raw_counts} ${meta.sample}_decontx.h5ad ${params.correct_ambient_RNA}
	"""
    
    stub:
    """
    touch ${meta.sample}_decontx.h5ad
    """	
}

 /*
 * create_h5ad
 */

 process basic_QC {
	tag "$meta.sample"
	container = "gcfntnu/scanpy:1.11.4"
	
	input:
	tuple val(meta), path(cellranger_h5), path(raw_counts)
	
	output:
	tuple val(meta), path("${meta.sample}.h5ad"), path(raw_counts)
	
	script:
    """
	export NUMBA_CACHE_DIR="./numbacache"
	create_h5ad.py ${cellranger_h5} ${meta.sample}.h5ad ${params.mitochondrial_prefix}
	"""
    
    stub:
    """
    touch ${meta.sample}.h5ad
    """	
}


/*
 * calculate doublet score
 */

 process scDblFinder {
	tag "$meta.sample"
	container = "tjmgibson/scrnaseq_preprocess:v2"
	publishDir "${params.results_dir}/h5ad_files/unfiltered/", mode: 'copy'
	
	input:
	tuple val(meta), path(h5ad)
	
	output:
	tuple val(meta), path("${meta.sample}_doublets.h5ad")
	
	script:
    """
	scDblFinder.R ${h5ad} ${meta.sample}_doublets.h5ad
	"""
    
    stub:
    """
    touch ${meta.sample}_doublets.h5ad
    """	
}

/*
 * generate QC plots
 */

 process QC_plots {
	tag "$meta.sample"
	container = "tjmgibson/scrnaseq_preprocess:v2"
	publishDir "${params.results_dir}/QC/preclustering/",pattern: '*.pdf', mode: 'copy'
	input:
	tuple val(meta), path(h5ad)
	path(util_file)
	
	output:
	tuple val(meta), path(h5ad), emit: adata
	path("${meta.sample}_precluster_QC.pdf"), emit: QC_plot

	
	script:
    """
	generate_QC_plots.R \
		${h5ad} \
		${meta.sample}_precluster_QC.pdf \
		${params.min_UMI} \
		${params.max_UMI} \
		${params.min_genes_per_cell} \
		${params.max_percent_mito} \
		${params.max_doublet_score} \
		${params.MAD_thresh} \
		${params.percentile_thresh} \
		${util_file}
		
	"""
    
    stub:
    """
    touch ${meta.sample}_precluster_QC.pdf
    """	
}

/*
 * perform QC prefiltering
 */

 process QC_prefilter {
	tag "$meta.sample"
	container = "tjmgibson/scrnaseq_preprocess:v2"
	publishDir "${params.results_dir}/h5ad_files/filtered/",pattern: '*.h5ad', mode: 'copy'
	publishDir "${params.results_dir}/filtering_stats/",pattern: '*.csv', mode: 'copy'
	
	input:
	tuple val(meta), path(h5ad)
	path(util_file)

	output:
	path("${meta.sample}_filtered.h5ad"), emit: filtered_adata
	path("${meta.sample}_precluster_QC_stats.csv"), emit: QC_stats

	
	script:
    """
	QC_prefiltering.R \
		${h5ad} \
		${meta.sample}_filtered.h5ad \
		${params.UMI_filter_strategy} \
		${params.gene_per_cell_filter_strategy} \
		${params.mito_filter_strategy} \
		${params.doublet_filter_strategy} \
		${params.min_UMI} \
		${params.max_UMI} \
		${params.min_genes_per_cell} \
		${params.max_percent_mito} \
		${params.max_doublet_score} \
		${params.MAD_thresh} \
		${params.percentile_thresh} \
		${params.R_utils_file} \
		${meta.sample}_precluster_QC_stats.csv
		
	"""
    
    stub:
    """
	touch ${meta.sample}_filtered.h5ad
    touch ${meta.sample}_precluster_QC.pdf
	touch ${meta.sample}_precluster_QC_stats.csv
    """	
}

process clustering {
	tag "${params.experiment_name}"
	container = "tjmgibson/scrnaseq_preprocess:v2"
	publishDir "${params.results_dir}/seurat_objects/", pattern: '*.rds', mode: 'copy'
	publishDir "${params.results_dir}/clusters/", pattern: '*clusters.csv', mode: 'copy'
	publishDir "${params.results_dir}/UMAP/", pattern: '*UMAP.csv', mode: 'copy'
	
	
	input:
	path(h5ad)
	
	output:
	path("${params.experiment_name}_clustered_seurat.rds"), emit: clustered_seurat
	path("*clusters.csv")
	path("*UMAP.csv")

	
	script:
	String h5ad_csv = h5ad.join(",")
    """
	clustering.R \
		${h5ad_csv} \
		${params.experiment_name}_clustered_seurat.rds \
		${params.cluster_npcs} \
		${params.cluster_resolution} \
		${params.integrate_samples} \
		${params.experiment_name} 
	"""
    
    stub:
    """
    touch ${params.experiment_name}_clustered_seurat.rds
	touch unintegrated_clusters.csv
	touch RPCA_integrated_clusters.csv
	touch unintegrated_UMAP.csv
	touch RPCA_UMAP.csv
    """	
}

/*
 * Run workflow
 */
 
workflow {
	
	cellranger_h5_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
	| splitCsv( header:true )
    | map { row ->
        meta = row.subMap('sample')
        [
        	meta, 
        	file(row.cellranger_h5, checkIfExists: true),
			file(row.cellranger_h5_raw, checkIfExists: true)]
    }
	
	h5ad_ch = decontX(cellranger_h5_ch)

	plot_ch = h5ad_ch
	| basic_QC
	| scDblFinder
	
	QC_ch = QC_plots(
		plot_ch,
		params.R_utils_file
	)

	filtered_adata_ch = QC_prefilter(
		QC_ch.adata,
		params.R_utils_file
		)

	filtered_adata_ch.filtered_adata
	| collect
	| clustering

}