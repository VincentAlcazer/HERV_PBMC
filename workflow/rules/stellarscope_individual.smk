rule cell_sort:
    conda:
        "../envs/samtools.yaml"
    output:
        "results/starsolo_alignment/{s}/{s}.Aligned.sortedByCB.bam"
    input:
        bam = "results/starsolo_alignment/{s}/{s}.collated.out.bam",
	barcodes = "results/starsolo_alignment/{s}/{s}.Solo.out/Gene/filtered/barcodes.tsv"
    log:
        "results/telescope_individual/{s}/cellsort.log"
    threads:
        config['telescope_threads']
    params:
        tmpdir = config['cell_sort_tmp']
    shell:
        '''
	samtools view -u -F 4 -D CB:{input.barcodes} {input.bam} | samtools sort -@ {threads} -n -t CB > {output}
	'''


rule stellarscope_individual:
    conda:
        "../envs/telescope.yaml"
    output:
        "results/telescope_individual/{s}/{s}_individual-TE_counts.mtx",
	directory("results/telescope_individual/{s}")
    input:
        bam = "results/starsolo_alignment/{s}/{s}.Aligned.sortedByCB.bam",
	annotation = rules.telescope_annotation.output,
	barcodes = "results/starsolo_alignment/{s}/{s}.Solo.out/Gene/filtered/barcodes.tsv"
    benchmark:
        "benchmarks/telescope_individual/{s}_telescope_individual.tsv"
    log:
        "results/telescope_individual/{s}/telescope.log"
    threads:
        config['telescope_threads']
    params:
        tmpdir = config['local_tmp'],
	out = "results/telescope_individual/{s}",
	exp_tag = "{s}_individual"
    shell:
        '''
	stellarscope assign\
	    --pooling_mode individual\
	    --updated_sam\
	    --exp_tag {params.exp_tag}\
	    --outdir {params.out}\
	    {input[0]}\
	    {input[1]}\
	    --barcodefile {input[2]}\
	    2>&1 | tee {log[0]}
	'''

