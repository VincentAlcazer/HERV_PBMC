rule index_genome_star:
    conda:
        "../envs/star.yaml"
    input:
        genome_fasta_gz = config['sequences']['genome_gz'],
	annotation_gtf_gz = config['annotations']['gencode_gz']
    output:
        directory(config['genome_index']['star'])
    params:
        sjdbOverhang = config['splice_junction_overhang']
    threads:
        config['star_index_threads']
    resources:
        mem_mb=config['star_index_mem_mb']
    shell:
        """
	# decompress fasta and gtf to a temporal directory
	tdir=$(mktemp -d {config[local_tmp]}/{rule}.XXXXXX)
	
	pigz -dc {input.genome_fasta_gz} > $tdir/genome.fa
	pigz -dc {input.annotation_gtf_gz} > $tdir/annotation.gtf
	
	# index reference genome with Star
	STAR\
	    --runThreadN {threads}\
	    --runMode genomeGenerate\
	    --genomeDir {output}\
	    --outFileNamePrefix {output}\
	    --genomeFastaFiles $tdir/genome.fa\
	    --sjdbGTFfile $tdir/annotation.gtf\
	    --sjdbOverhang {params.sjdbOverhang}
	"""

