rule starsolo_alignment:
    """
    Align sequencing reads from 10x V3 single-cell RNA-seq experiments
    """
    conda:
        "../envs/star.yaml"
    input:
        fastq_files = "data/{s}_fastqs",
	genome = config['genome_index']['star'],
	whitelist = config['whitelist']['v3']
    output:
        "results/starsolo_alignment/{s}/{s}.Aligned.sortedByCoord.out.bam"
    params:
        cDNA_reads = lambda wc: ','.join(samples.loc[wc.s]["read2"]), # reverse reads (10x v3)
	CB_reads = lambda wc: ','.join(samples.loc[wc.s]["read1"]), # forward reads (10x v3)
	out_prefix = "results/starsolo_alignment/{s}/{s}.",
	cb_start = config['cellbarcode_start'],
	cb_length = config['cellbarcode_length'],
	umi_start = config['umi_start'],
	umi_length = config['umi_length'],
	max_multimap = config['max_multimap']
    benchmark:
        "benchmarks/starsolo_alignment/{s}_starsolo_alignment.tsv"
    threads:
        config['star_alignment_threads']
    shell:
        '''
	#--- STARsolo (turned on by --soloType CB_UMI_Simple)
	STAR\
	    --runThreadN {threads}\
	    --genomeDir {input.genome}\
	    --readFilesIn {params.cDNA_reads} {params.CB_reads}\
	    --readFilesCommand gunzip -c\
	    --soloType CB_UMI_Simple\
	    --soloCBwhitelist {input.whitelist}\
	    --soloCBstart {params.cb_start}\
	    --soloCBlen {params.cb_length}\
	    --soloUMIstart {params.umi_start}\
	    --soloUMIlen {params.umi_length}\
	    --outFilterMultimapNmax {params.max_multimap}\
	    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM\
	    --outSAMtype BAM SortedByCoordinate\
	    --outFileNamePrefix {params.out_prefix}
	'''

