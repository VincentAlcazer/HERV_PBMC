checkpoint split_bam_in_clusters:
    conda:
        "../envs/telescope.yaml"
    output:
        touch("results/starsolo_alignment/{s}/.{s}_split_bam")
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

