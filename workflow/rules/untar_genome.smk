# I don't love this rule
# maybe we should keep the fasta as tar gz 
# and simply decompress it to temp 
# in the rule that indexes the genome
#
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
