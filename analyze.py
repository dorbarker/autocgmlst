#!/usr/bin/env python3

from math import log
import pandas as pd
import subprocess

def shannon(sequence, base=2):
   
    def p(i):
        return sequence.count(i) / len(sequence)

    return -sum(p(i) * log(p(i), base) for i in set(sequence))

def shannon_by_call():
    '''Calculates Shannon Entropy for the relative abundances
    of allele calls for each locus
    '''
    pass

def shannon_by_sequence():
    '''Calculates Shannon Entropy for a particular k-mer of an
    aligned sequence
    '''
    pass

def map_contig_locations():
    pass

def adjusted_wallace():
    pass

def precompute_clusters():
    pass

def subset_genes():
    pass

