rule split_bam_clusters:
    conda:
        "../envs/telescope.yaml"
    output:
        "results/starsolo_alignment/{s}/{s}.Aligned.sortedByCB.CLUST_1.bam"
    input:
        bam = "results/starsolo_alignment/{s}/{s}.Aligned.sortedByCB.bam",
        clusters = "data/{s}_clusters/analysis/clustering/graphclust/clusters.csv"
    benchmark:
        "benchmarks/split_bam_clusters/{s}_split_clusters.tsv"
    threads:
        config['samtools_collate_threads']
    shell:
        '''
        python workflow/scripts/split_by_cluster.py {input.bam} {input.clusters}
        '''


rule obtain_cluster_barcodes:
    """
    Notice that we have to use double curly braces for awk and escape the quotes in gsub 
    """
    output:
        "results/starsolo_alignment/{s}/{s}.Aligned.sortedByCB.CLUST_1.barcodes"
    input:
        "data/{s}_clusters/analysis/clustering/graphclust/clusters.csv"
    params:
        filename="results/starsolo_alignment/{s}/{s}.Aligned.sortedByCB"
    shell:
        """
        file={params.filename}
        # awk receives a file (the rule input), the variable f (the rule params) 
        # awk skips the first record, substitutes the -1 suffix with nothing in the first column
        # prints the formatted first column to the file whose name is built using the second column
        #
        awk -F, -v f="$file" 'NR>1 {{gsub(\"-1\",\"\",$1); print $1 > f".CLUST_"$2".barcodes"}}' {input}
        """

