#! /usr/bin/env python
# -*- coding utf-8 -*-

################################# DOWNLOAD REFS #################################

rule download_remote:
    """ Downloads a remote file and checks the md5sum
    """
    output:
        'resources/downloads/{file}'
    params:
        url = lambda wildcards: config['downloads'][wildcards.file]['url'],
        md5 = lambda wildcards: config['downloads'][wildcards.file]['md5']
    shell:
        '''
	curl -L {params.url} > {output}
	echo {params.md5} {output} | md5sum -c -
        '''


################################# EXTRACT REFS #################################
#
# this is a highly specific rule
#rule untar_genome:
#    input:
#        'resources/downloads/GRCh38.d1.vd1.fa.tar.gz'
#    output:
#        config['sequences']['genome_gz']#,
#        config['sequences']['genome_idx'],
#        config['sequences']['genome_dict']
#    conda:
#        "../envs/utils.yaml"
#    shell:
#        '''
#	#mkdir -p $(dirname {output[0]})
#	tar -Oxzf {input} | bgzip > {output}
#	#samtools faidx {output[0]}
#	#picard CreateSequenceDictionary R={output[0]} O={output[2]}
#        '''
#
#
#rule extract_transcriptome:
#    input:
#        'resources/downloads/gencode.v38.annotation.gtf.gz',
#        config['sequences']['genome']
#    output:
#        config['sequences']['transcripts'],
#        config['sequences']['transcripts_dupinfo'],
#        config['sequences']['transcripts_list']
#    conda:
#        "../envs/utils.yaml"
#    shell:
#        '''
#tfa=$(mktemp -p {config[local_tmp]})
#gunzip -c {input[1]} > $tfa
#gunzip -c {input[0]} | gffread - -M -d {output[1]} -g $tfa -w {output[0]}
#grep ">" {output[0]} | sed "s/^>//" | sort | uniq > {output[2]}
#rm -f $tfa*
#        '''

################################## ID MAPPING ##################################

#rule id_mapping:
#    input:
#        'refs/downloads/gencode.v38.annotation.gtf.gz',
#        config['sequences']['transcripts_list']
#    output:
#        config['annotations']['ttg'],
#        config['annotations']['gsym']
#    run:
#        tx_g = {}
#        g_sym = {}
#
#        raw = (bl.decode('utf-8') for bl in gzip.open(input[0], 'rb'))
#        lines = (r.strip('\n').split('\t') for r in raw if not r.startswith('#'))
#        for l in lines:
#            d = dict(t for t in re.findall('(\S+)\s+"([\s\S]*?)";', l[8]))
#            if l[2] == 'gene':
#                if d['gene_id'] in g_sym:
#                    assert g_sym[d['gene_id']] == d['gene_name'], "Gene name mismatch: %s %s" % (d['gene_name'], g_sym[d['gene_id']])
#                g_sym[d['gene_id']] = d['gene_name']
#
#            if l[2] == 'transcript':
#                if d['transcript_id'] in tx_g:
#                    assert tx_g[d['transcript_id']] == d['gene_id'], "Gene ID mismatch: %s %s" % (d['gene_id'], tx_g[d['transcript_id']])
#                tx_g[d['transcript_id']] = d['gene_id']
#
#        # Generate ttg
#        with open(output[0], 'w') as outh:
#            print('TXNAME\tGENEID', file=outh)
#            txlist = (l.strip('\n').split()[0] for l in open(input[1], 'rU'))
#            for tx in txlist:
#                print('%s\t%s' % (tx, tx_g[tx]), file=outh)
#
#        # Generate gsym
#        with open(output[1], 'w') as outh:
#            print('GENEID\tSYM', file=outh)
#            for t in g_sym.items():
#                print('%s\t%s' % t, file=outh)
#
############################### FORMAT ANNOTATIONS ##############################
#
#rule telescope_annotation:
#    input:
#        'refs/downloads/retro.hg38.v1.gtf',
#        config['sequences']['genome_idx']
#    output:
#        config['annotations']['retro']
#    shell:
#        '''
#python workflow/scripts/sortgtf.py --fai {input[1]} < {input[0]} > {output[0]}
#        '''
#
#
################################# RULE TARGETS #################################
#
#rule complete_download:
#    input:
#        config['sequences']['genome'],
#        config['sequences']['genome_idx'],
#        config['sequences']['genome_dict'],
#        config['sequences']['transcripts'],
#        config['sequences']['transcripts_dupinfo'],
#        config['sequences']['transcripts_list'],
#        config['annotations']['ttg'],
#        config['annotations']['gsym']

