import argparse
import os
import json
import pandas as pd
from utilities import Sequence
from collections import Counter, defaultdict
from glob import glob
from Bio import SeqIO
from itertools import combinations, product
import operator

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--alleles', required=True, metavar='PATH',
                        help='Directory containing FASTA-formatted \
                              cgMLST alleles')

    parser.add_argument('--calls', required=True, metavar='PATH',
                        help='Path to table of allele calls')

    parser.add_argument('--jsons', required=True, metavar='PATH',
                        help='Directory containing MIST output in \
                              JSON format')

    return parser.parse_args()

def load_alignments(gene_dir):

    genes = {}

    for g in glob(gene_dir + '/*'):

        bn, _ = os.path.splitext(os.path.basename(g))

        with open(g, 'r') as f:
            alleles = {}
            for rec in SeqIO.parse(f, 'fasta'):

                s = str(rec.seq).replace('-', '')
                #alleles[rec.id] = Sequence(s)
                alleles[rec.id] = s
            genes[bn] = alleles

    return genes

def allele_abundances(locus, calls):

    return Counter((x for x in calls[locus] if x >= 1))

def match_partial_alleles(fragment, alleles):

    return {k: v for k, v in alleles.items() if fragment in v}

def partial_match_abundances(abundances, partial_matches):

    return {k: abundances[k] for k in partial_matches}

def combine_probabilities(links, partial_abundances):


def linkage(centre, calls):

    associations = defaultdict(dict)

    ind = calls.columns.get_loc(centre)

    # wraps around the table
    left = calls.columns[ind - 1]
    right = calls.columns[ind + 1 if ind + 1 < len(calls.columns) else 0]


    triplet = calls[[left, centre, right]]
    #condition = triplet[centre] >= 1 & triplet[left] >= 1 & triplet[right] >= 1
    condition = [all(map(lambda x: x >= 1, t)) for t in triplet.values]
    clean = triplet[condition]
    #flanks = clean[[left, right]]

    for l, c, r in clean.values:

        try:
            associations[ (l, r) ][c] += 1
        except KeyError:
            associations[ (l, r) ][c] = 1

    return associations

def recover_gene_fragment(gene, genome, jsons):

    with open(os.path.join(jsons, genome + '.json'), 'r') as f:
        data = json.load(f)

        results = data['Results'][0]['TestResults']

        for key in results:
            if gene in results[key]:
                try:
                    s = results[key][gene]['BlastResults']['SubjAln']
                #return Sequence(s)
                    return s
                except TypeError:
                    return 'X'

def predict_allele(locus, genome, genes, jsons, calls):

    alleles = genes[locus]

    fragment = gene_fragment(locus, genome, jsons)

    matches = match_partial_alleles(fragment, alleles)

    associations = linkage(locus, calls)

def rarefy(n, locus, calls):

    def n_comb(N, n):
        return len(combinations(N, n))
    non_missing = list(filter(lambda x: x >= 1, calls[locus]))

    N = len(non_missing)
    K = set(non_missing)

    return len(K) - (n_comb(N, n)**-1) * sum(n_comb(N - non_missing.count(i), n) for i in K)

def calculate_probabilities(gene_name, gene_alleles, genome, abundances,
                            neighbour_calls, associations, jsons):
    try:
        links = associations[neighbour_calls]

    except KeyError:
        links = {}

    fragment = recover_gene_fragment(gene_name, genome, jsons)

    partial_matches = match_partial_alleles(fragment, gene_alleles)

    partial_abundances = partial_match_abundances(abundances, partial_matches)

    fragment_probabilities = calculate_fragment_probabilities(partial_abundances,
                                                              novel_allele_probability)

    flanking_probabilities = calculate_flanking_probabilities(partial_matches, links)

    combined_probabilities = combine_probabilities(links, partial_abundances)

def calculate_fragment_probabilities(fragment_abundances, novel_allele_probability):

    total = sum(fragment_abundances.values())
    l = len(fragment_abundances)

    if total is 0:
        # subtract novel_allele_probability
        return {k: 1 / l for k in fragment_abundances}

    elif any(x is 0 for x in fragment_abundances.values()):

        # calculate probabilities for represented alleles
        # subtract novel_allele_probability equally from each
        # set any zeros to novel_allele_probability


        pass

    else:
        return {k: v / total for k, v in fragment_abundances.items()}


def calculate_flanking_probabilities(partial_matches, links):

    possible_links = {k: v for k, v in links.items() if k in partial_matches}

    total = sum(possible_links.values())
    l = len(possible_links)


    return {k: v / total for k, v in possible_links.items()}

def fix_allele_calls(calls, alleles, jsons):

    out_calls = calls.copy(deep=True)

    for gene in calls:

        position = calls.columns.get_loc(gene)

        neighbours = [calls.columns[position - 1], calls.columns[position + 1]]

        affected_genomes = calls.index[calls[gene] < 1]

        abundances = allele_abundances(gene, calls)

        associations = linkage(gene, calls)

        gene_alleles = alleles[gene]

        for genome in affected_genomes:

            neighbour_calls = tuple(calls[neighbours].loc[[genome]])

            probabilities = calculate_probabilities(gene, gene_alleles, genome, abundances,
                                                    neighbour_calls, associations, jsons)

            #out_calls[gene][genome] = best

    return out_calls

def main():

    args = arguments()

    genes = load_alignments(args.alleles)

    calls = pd.read_csv(args.calls, index_col='genomes')

    fixed = fix_allele_calls(calls, genes, args.jsons)


if __name__ == '__main__':

    main()
