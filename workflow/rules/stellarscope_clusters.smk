#! /usr/bin/env python
# -*- coding utf-8 -*-
import os

# function to get the clusters bams
def get_bam_files(wildcards):
    bams_dir = os.path.dirname(checkpoints.split_bam_in_clusters.get(**wildcards).output[0])
    SAMPLE, CLUSTER = glob_wildcards(os.path.join(bams_dir, "{S}.Aligned.sortedByCB.CLUST_{C}.bam"))
    return expand(os.path.join(bams_dir, "{{s}}.Aligned.sortedByCB.CLUST_{c}.bam"), c=CLUSTER)

# function to get the clusters barcodes
def get_barcode_files(wildcards):
    bc_dir = os.path.dirname(checkpoints.split_barcodes_in_clusters.get(**wildcards).output[0])
    SAMPLE, CLUSTER = glob_wildcards(os.path.join(bc_dir, "{S}.CLUST_{C}.barcodes"))
    return expand(os.path.join(bc_dir, "{{s}}.CLUST_{c}.barcodes"), c=CLUSTER)

rule stellarscope_clusters:
    conda:
        "../envs/telescope.yaml"
    output:
        "results/telescope_clusters/{s}/{s}_clusters_CLUST_{c}-TE_counts.mtx"
    input:
        get_bam_files,
        get_barcode_files, 
	annotation = rules.telescope_annotation.output
    benchmark:
        "benchmarks/telescope_clusters/{s}_telescope_clusters_CLUST_{c}.tsv"
    log:
        "results/telescope_clusters/{s}/telescope_CLUST_{c}.log"
    threads:
        config['telescope_threads']
    params:
        bam = "results/starsolo_alignment/{s}/{s}.Aligned.sortedByCB.CLUST_{c}.bam",
        barcodes = "data/{s}_clusters/analysis/clustering/graphclust/{s}.CLUST_{c}.barcodes",
        tmpdir = config['local_tmp'],
	out = "results/telescope_clusters/{s}",
	exp_tag = "{s}_clusters_CLUST_{c}"
    shell:
        '''
	getrss() {{
	    local cmd=$1
	    local param=$2
	    ps -C $cmd -o args= | grep "$param" | awk '{{$1=int(100 * $1/1024/1024)/100"GB";}}{{ print $1;}}' | while read v; do echo "Memory usage (RSS) for $cmd (param: $param): $v"; done
	}}

	while true; do getrss stellarscope {wildcards.s}; sleep 5; done &

	stellarscope assign\
	    --updated_sam\
	    --exp_tag {params.exp_tag}\
	    --outdir {params.out}\
	    {params.bam}\
	    {input.annotation}\
	    --barcodefile {params.barcodes}\
	    2>&1 | tee {log[0]}
	'''

