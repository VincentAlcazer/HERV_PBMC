#!/usr/bin/env python
import pysam
from collections import Counter

#alignment_bam = snakemake.input["bam"]
#cluster_csv = snakemake.input["clusters"]

def split_cluster(alignment_bam, cluster_csv, verbose=False):
    """ Parse the cluster file """
    lines = (l.strip().split(',') for l in open(cluster_csv, 'r'))
    header = next(lines)
    mapping = {}
    clusters = set()
    for cb,cl in lines:
        # Cell barcodes from 10x analysis have '-1' suffix while STARsolo BAM files
        # have only the barcode sequence. Splitting the barcode on '-' should have 
        # no effect if CSV does not use the suffix.
        mapping[cb.split('-')[0]] = cl
        clusters.add(cl)
    
    print('%d clusters found' % len(clusters))
    print('%d cell barcodes found' % len(mapping))
    
    """ Create output BAM files """
    bamfile = pysam.AlignmentFile(alignment_bam, "rb")
    output_bams = {}
    for cl in sorted(clusters):
        fn = '%s.CLUST_%s.bam' % ('.'.join(alignment_bam.split('.')[:-1]), cl)
        print('creating output BAM for CLUST_%s: %s' % (cl, fn))
        output_bams[cl] = pysam.AlignmentFile(fn, "wb", template=bamfile)
    
    """ Fetch reads """
    n_alignments = 0
    n_fail = 0
    aln_per_clust = Counter()
    for r in bamfile.fetch(until_eof=True):
        n_alignments += 1
        if not n_alignments % 1000000:
            print('%d alignments processed' % n_alignments)
        try:
            aln_bc = r.get_tag('CB')
        except KeyError:
            if verbose: print('CB tag not found')
            n_fail += 1
            continue
        
        try:
            dest = mapping[aln_bc]
        except KeyError:
            if verbose: print('Barcode %s not in clusters' % aln_bc)
            n_fail += 1        
            continue
        
        if output_bams[dest].write(r):
            aln_per_clust[dest] += 1
    
    bamfile.close()
    for fh in output_bams.values():
        fh.close()
    
    print('%d total alignments' % n_alignments)
    print('%d alignments excluded' % n_fail)
    for cl in sorted(clusters):
        print('CLUST_%s: %d alignments' % (cl, aln_per_clust[cl]))
    
    return [(output_bams[cl].filename, aln_per_clust[cl]) for cl in sorted(clusters)]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Split BAM file by cluster.')
    parser.add_argument('alignment_bam',
                        help='BAM file to be split.'
                        )
    parser.add_argument('cluster_csv',
                        help='CSV file with cell barcode in column 1, cluster ID in column 2'
                        )
    args = parser.parse_args()
    rv = split_cluster(args.alignment_bam, args.cluster_csv)
