rule get_whitelist:
    """
    Get appropriate 10x whitelist of used barcodes
    """
    input:
        "resources/downloads/whitelist.10x.v3.txt.gz"
    output:
        config['whitelist']['v3']
    shell:
        '''
        gunzip -c {input} > {output}
        '''

