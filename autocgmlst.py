#!/usr/bin/env python3

from functools import partial
from multiprocessing import cpu_count
from utilities import contents, fix_headers
import argparse
import os
import subprocess
import scheme_create

def arguments():

    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()

    ### Create
    create = subparsers.add_parser('create', help='Create a cgMLST scheme')

    create.add_argument('--genomes', required=True, help='Genome directory')

    create.add_argument('--work-dir', required=True, help='Working directory')

    create.add_argument('--identity', type=float, required=True,
                        help='Minimum sequence identity \
                              for homologue determination')

    create.add_argument('--coverage', type=float, required=True,
                        help='Minimum sequence coverage for homology')

    create.add_argument('--cores', type=int, default=cpu_count(),
                        help='Number of CPU cores to use')

    create.add_argument('--mist', help='Path to MIST.exe')

    create.add_argument('--refine-identity', type=float,
                        help='Refine homologues with this identity threshold')

    create.add_argument('--refine-coverage', type=float,
                        help='Refine homologues with this coverage threshold')

    create.add_argument('--genome-quality-cutoff', type=float,
                        help='Maximum number of missing loci permitted in a genome')

    create.set_defaults(func=scheme_create.create)

    ### Analyze
    analyze = subparsers.add_parser('analyze', help='Analyze a cgMLST scheme')

    return parser.parse_args()

def prepare_create_args(args):

    subdir = partial(os.path.join, args.work_dir)

    fargs = {'genome_dir': subdir('genomes'),
             'alleles_dir': subdir('alleles'),
             'json_dir': subdir('jsons'),
             'work_dir': args.work_dir,
             'prokka_out': subdir('prokka_out'),
             'min_identity': args.identity,
             'min_coverage': args.coverage,
             'refine_identity': args.refine_identity or args.identity,
             'refine_coverage': args.refine_coverage or args.coverage,
             'mist_bin': args.mist,
             'genome_quality_cutoff': args.genome_quality_cutoff,
             'cores': args.cores}

    for i in ('work_dir', 'genome_dir', 'alleles_dir', 'json_dir', 'prokka_out'):
        if not os.access(fargs[i], os.F_OK):
            os.mkdir(fargs[i])

    for j in contents(args.genomes):
        o = os.path.join(fargs['genome_dir'], os.path.basename(j))
        fix_headers(j, o)

    return fargs

def main():

    args = arguments()

    if args.func is scheme_create.create:

        fargs = prepare_create_args(args)

    args.func(**fargs)

if __name__ == '__main__':
    main()
