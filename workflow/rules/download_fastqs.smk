rule download_fastqs:
    """ Download fastq files from 10x
    """
    output:
        "samples/{s}/{s}_fastqs.tar"
    benchmark:
        "benchmarks/download_fastqs/{s}_download_fastqs.tsv"
    params:
        url = lambda wc: samples[wc.s]['url'],
	md5 = lambda wc: samples[wc.s]['md5']
    shell:
        '''
	curl -L {params.url} > {output}
	echo {params.md5} {output} | md5sum -c -

	'''

rule untar_fastqs:
    input:
        "samples/{s}/{s}_fastqs.tar"
    output:
        "samples/{s}/{s}_fastqs"
    shell:
        '''
	tar -xvf {input}

rule fastqs_complete:
    input:
        expand("samples/{s}/{s}_fastqs", s = samples[sample_name])

