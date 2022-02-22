rule download_clusters:
    """ Download graphclust files from 10x
    """
    output:
        "data/{s}_analysis.tar.gz"
    benchmark:
        "benchmarks/download_clusters/{s}_download_clusters.tsv"
    params:
        url = lambda wc: samples.loc[wc.s]["url_clusters"],
	md5 = lambda wc: samples.loc[wc.s]["md5_clusters"]
    shell:
        '''
	curl -L {params.url} > {output}
	echo {params.md5} {output} | md5sum -c -

	'''

rule untar_clusters:
    input:
        "data/{s}_analysis.tar.gz"
    output:
        "data/{s}_clusters/analysis/clustering/graphclust/clusters.csv"
    params:
        outdir = "data/{s}_clusters"
    shell:
        '''
	tar -xvf {input} --directory {params.outdir}
	'''

#rule fastqs_complete:
#    input:
#        expand("data/{s}/{s}_fastqs", s = samples[sample_name])

