rule download_fastqs:
    """ Download fastq files from 10x
    """
    output:
        "data/{s}/{s}_fastqs.tar"
    benchmark:
        "benchmarks/download_fastqs/{s}_download_fastqs.tsv"
    params:
        url = lambda wc: samples.loc[wc.s]["url"],
	md5 = lambda wc: samples.loc[wc.s]["md5"]
    shell:
        '''
	curl -L {params.url} > {output}
	echo {params.md5} {output} | md5sum -c -

	'''

rule untar_fastqs:
    input:
        "data/{s}/{s}_fastqs.tar"
    output:
        "data/{s}/{s}_fastqs"
    params:
        outdir = "data/{s}"
    shell:
        '''
	tar -xvf {input} --directory {params.outdir}
	'''

#rule fastqs_complete:
#    input:
#        expand("data/{s}/{s}_fastqs", s = samples[sample_name])

