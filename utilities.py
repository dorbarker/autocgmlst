from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from functools import partial
import collections
import glob
import operator
import os
import subprocess

Gene = collections.namedtuple('Gene', ['id', 'seq'])

def basename(path):

    name, _ = os.path.splitext(os.path.basename(path))
    return name

def calculate_distance(s1, s2):

    def aln_parse(aln):

        lines = aln.strip().split()

        sec_head = lines.index('>seq2')

        return ''.join(lines[1:sec_head]), ''.join(lines[sec_head + 1:])

    def format_input(seq1, seq2):
        return '>seq1\n{}\n>seq2\n{}\n'.format(seq1, seq2)

    def hamming(a, b):

        shared = [x for x in zip(a, b) if '-' not in x]
        len_shared = len(shared)

        identity = 1 - (sum(int(i != j) for i, j in shared) / len_shared)
        cov = len(a) / len_shared

        return identity, cov

    mafft = ('mafft', '--quiet', '--retree', '1', '-')
    result = subprocess.check_output(mafft, input=format_input(s1, s2),
                                     universal_newlines=True)

    return hamming(*aln_parse(result))

def contents(path, pattern=''):

    return glob.iglob(os.path.join(path, '*' + pattern))

def load_gene(record):

    return record.id, str(record.seq)

def fix_headers(infile, outfile):
    
    def fix_header(i, record):
        o = record
        o.id = str(i)
        return o
    
    with open(infile, 'r') as i, open(outfile, 'w') as o:
        records = (fix_header(a, b) for a, b in enumerate(SeqIO.parse(i, 'fasta'), 1))
        
        SeqIO.write(records, o, 'fasta')

def match_gene(g2, g1, min_identity, min_coverage, gene_dict):

    ident, cov = calculate_distance(gene_dict[g1],
                                    gene_dict[g2])

    if ident > min_identity and cov > min_coverage:
        return g2

def match_genes(gene_dict, min_identity, min_coverage, cores):
    
    def to_search(g2, g1, lengths, to_skip):

        return  g2 not in to_skip and lengths[g2] / lengths[g1] >= min_coverage

    from time import time

    homologues = collections.defaultdict(list)

    lengths = {k: len(v) for k, v in gene_dict.items()}
    print(len(lengths))
    to_skip = set()
    t1 = time()
    starttime = t1
    # run the shortest genes first to eliminate as many search possibilities as possible
    # for the later, longer comparisons
    counter = 0
    
    for gene1 in sorted(gene_dict, key=lambda x: lengths[x], reverse=True):
        counter += 1
        
        if gene1 in to_skip:
            continue

        with ThreadPoolExecutor(cores) as tpe:

            f = partial(match_gene, g1=gene1, min_identity=min_identity, min_coverage=min_coverage, gene_dict=gene_dict)
            
            filter_func = partial(to_search, g1=gene1, lengths=lengths, to_skip=to_skip | {gene1})
            
            l = filter(filter_func, gene_dict)
            
            results = tpe.map(f, l)

            for r in filter(None, results):
                homologues[gene1].append( [r] + homologues[r] )
                del homologues[r]
                to_skip.add(r)

        t2 = time()
        
        print(counter, len(homologues), len(to_skip), t2-t1, t2-starttime)
        t1 = t2
    return homologues
