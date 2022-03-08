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

checkpoint split_barcodes_in_clusters:
    """
    Notice that we have to use double curly braces for awk and escape the quotes in gsub
    """
    output:
        touch("data/{s}_clusters/analysis/clustering/graphclust/.{s}_split_barcodes")
    input:
        "data/{s}_clusters/analysis/clustering/graphclust/clusters.csv"
    params:
        name="data/{s}_clusters/analysis/clustering/graphclust/{s}"
    shell:
        """
        file={params.name}
        # awk receives a csv file (the rule input), the variable f (the rule params)
        # awk skips the first record, substitutes the -1 suffix with nothing in the first column
        # then prints this formatted first column to a file whose name is built using the second column
        #
        awk -F, -v f="$file" 'NR>1 {{gsub(\"-1\",\"\",$1); print $1 > f".CLUST_"$2".barcodes"}}' {input}
	"""

