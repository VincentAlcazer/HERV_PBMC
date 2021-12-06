#! /usr/bin/env python
# -*- coding utf-8 -*-

rule stellarscope:
    conda:
        "../envs/telescope.yaml"
    output:
        "results/telescope/{s}/{s}-TE_counts.tsv",
	"results/telescope/{s}/{s}-updated.bam",
	"results/telescope/{s}/{s}-other.bam"
    input:
        bam = "results/starsolo_alignment/{s}/{s}.collated.out.bam",
	annotation = rules.telescope_annotation.output
    benchmark:
        "benchmarks/telescope/{s}_telescope.tsv"
    log:
        "results/telescope/{s}/telescope.log"
    threads:
        config['telescope_threads']
    params:
        tmpdir = config['local_tmp'],
	out = "results/telescope/{s}",
	exp_tag = "{s}"
    shell:
        '''
	telescope sc assign\
	    --updated_sam\
	    --exp_tag {params.exp_tag}\
	    --outdir {params.out}\
	    {input[0]}\
	    {input[1]}\
	    2>&1 | tee {log[0]}
	'''

rule sample_complete:
    input:
        rules.stellarscope.output
    output:
        touch("results/completed/{s}_completed.txt")

