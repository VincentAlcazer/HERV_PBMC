#! /usr/bin/env python
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
import argparse

def sortgtf(input, chrom_order=None):
    gtf = [l.strip('\n').split('\t') for l in input if not l.startswith('#')]
    if chrom_order is None:
        gtf.sort(key=lambda x:int(x[3]))
        gtf.sort(key=lambda x:x[0])
    else:
        gtf = [g for g in gtf if g[0] in chrom_order]
        gtf.sort(key=lambda x:int(x[3]))
        gtf.sort(key=lambda x:chrom_order[x[0]])

    return gtf

def main(args):
    if args.sdict is None and args.fai is None:
        gtf = sortgtf(args.input)
    elif args.sdict:
        chroms = []
        dfile = (l.strip('\n').split('\t') for l in args.sdict)
        for r in dfile:
            if r[0] == '@SQ':
                d = {_[:2]:_[3:] for _ in r[1:]}
                chroms.append(d['SN'])
        chrom_order = {v:i for i,v in enumerate(chroms)}
        gtf = sortgtf(args.input, chrom_order)
    elif args.fai:
        chroms = [l.strip('\n').split('\t')[0] for l in args.fai]
        chrom_order = {v:i for i,v in enumerate(chroms)}
        gtf = sortgtf(args.input, chrom_order)

    print('\n'.join('\t'.join(_) for _ in gtf), file=args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sort GTF.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--sdict',
                        type=argparse.FileType('r'),
                        help='''Sequence dictionary (i.e. from picard
                                CreateSequenceDictionary) for chromosome order'''
    )
    group.add_argument('--fai',
                        type=argparse.FileType('r'),
                        help='''FASTA index with chromosome name in first column (i.e.
                                from samtools faidx) for chromosome order'''
    )
    parser.add_argument('input',
                        type=argparse.FileType('r'),
                        nargs='?',
                        default=sys.stdin,
                        help='GTF file to be sorted'
    )
    parser.add_argument('output',
                        type=argparse.FileType('w'),
                        nargs='?',
                        default=sys.stdout,
                        help='Sorted GTF file'
    )
    try:
        main(parser.parse_args())
    except BrokenPipeError:
        pass
