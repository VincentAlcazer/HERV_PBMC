rule samtools_collate:
    conda:
        "../envs/samtools.yaml"
    output:
        "results/starsolo_alignment/{s}/{s}.collated.out.bam"
    input:
        "results/starsolo_alignment/{s}/{s}.Aligned.sortedByCoord.out.bam"
    benchmark:
        "benchmarks/samtools_collate/{s}_samtools_collate.tsv"
    params:
        tmpdir = config['tmpdir'] 
    threads:
        config['samtools_collate_threads']
    shell:
        '''
	samtools collate\
	    {input}\
	    -o {output}\
	    -@ {threads}\
	    {params.tmpdir}
	'''

