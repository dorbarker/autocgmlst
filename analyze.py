#!/usr/bin/env python3

from math import log
import pandas as pd
import subprocess

def shannon(sequence, base=2):
   
    def p(i):
        return sequence.count(i) / len(sequence)

    return -sum(p(i) * log(p(i), base) for i in set(sequence))

def shannon_by_call(gene_calls):
    '''Calculates Shannon Entropy for the relative abundances
    of allele calls for each locus
    '''
    
    intact = gene_calls.apply(intact_alleles)
    
    shannons = map(shannon, intact)
    
    return shannons 

def shannon_by_sequence():
    '''Calculates Shannon Entropy for a particular k-mer of an
    aligned sequence
    '''
    pass

def intact_alleles(gene):

    return tuple(filter(lambda x: x not in (-1, 0), gene))

def map_contig_locations():
    pass

def adjusted_wallace():
    pass

def precompute_clusters():
    pass

def subset_genes(n, reps):
    pass

def failure_cost(gene):
    
    present = intact_alleles(gene) 
    
    N = len(gene)

    missing_probability = (N - len(present)) / N)

    H = shannon(present)

    return H - missing_probability * H
