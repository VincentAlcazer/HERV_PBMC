rule telescope_annotation:
    input:
        "resources/downloads/retro.hg38.v1.gtf",
	config['sequences']['genome_idx']
    output:
        config['annotations']['retro']
    shell:
        '''
	python workflow/scripts/sortgtf.py --fai {input[1]} < {input[0]} > {output[0]}
	'''

