rule untar_genome:
    conda:
        "../envs/decompress.yaml"
    input:
        config['sequences']['genome_tar_gz']
    output:
        config['sequences']['genome_gz']
    shell:
        """
        tar --extract --to-stdout --gzip --file {input} | pigz > {output}
        """
