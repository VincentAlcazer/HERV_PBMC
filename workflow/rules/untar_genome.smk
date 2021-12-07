rule untar_genome:
    conda:
        "../envs/decompress.yaml"
    input:
        config['sequences']['genome_tar_gz']
    output:
        config['sequences']['genome_gz'],
	config['sequences']['genome_idx']
    shell:
        """
        tar --extract --to-stdout --gzip --file {input} | bgzip > {output[0]}
	samtools faidx {output[0]}
        """
