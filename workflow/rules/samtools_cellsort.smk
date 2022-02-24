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

