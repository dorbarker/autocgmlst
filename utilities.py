from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import collections
import glob
import operator
import os
import subprocess
import tempfile
import re

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

def fasta_formatter(gene_dict, *names):

    seqs = ('>{}\n{}'.format(name, gene_dict[name]) for name in names)

    return '\n'.join(seqs)

def parse_cdhit_clusters(path, lengths):

    pt = re.compile('>.*\.\.\.')

    matching_clusters = collections.defaultdict(list)

    with open(path, 'r') as f:
        data = f.read().split('>Cluster')

    clusters = ((x.strip('>.') for x in re.findall(pt, y)) for y in data if y)

    size_sorted_clusters = (sorted(c, key=lambda n: -lengths[n]) for c in clusters)

    return {z[0]: set(z[1:]) for z in size_sorted_clusters}


def cluster_collapse(fasta, lengths, min_identity, min_coverage):

    with tempfile.TemporaryDirectory() as d:

        fasta_path = os.path.join(d, 'clusters.fasta')
        out_path = os.path.join(d, 'clusters.out')

        with open(fasta_path, 'w') as f:
            f.write(fasta)

        try:
            cdhit = 'cd-hit'

            subprocess.check_call(('which', cdhit), stdout=subprocess.DEVNULL)

        except CalledProcessError:

            cdhit = 'cdhit'

        cmd = (cdhit, '-i', fasta_path, '-o', out_path, '-c', str(min_identity),
               '-s', str(min_coverage), '-T', '1', '-M', '0', '-d', '0')

        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        return parse_cdhit_clusters(out_path + '.clstr', lengths)

def to_search(g2, g1, lengths, to_skip, min_identity, min_coverage):

    return  g2 not in to_skip and 1.0 >= (lengths[g2] / lengths[g1]) >= min_coverage

def match_genes(gene_dict, min_identity, min_coverage, cores):
    from time import time

    to_skip = set()

    lengths = {k: len(v) for k, v in gene_dict.items()}
    genes = sorted(gene_dict, key=lambda x: lengths[x], reverse=True)

    segments = [genes[i:i + cores] for i in range(0, len(genes), cores)]
    start = time()
    with ProcessPoolExecutor(cores) as ppe:
        f = partial(collapse_segment, gene_dict=gene_dict, lengths=lengths,
                    min_identity=min_identity, min_coverage=min_coverage)

        segment_results = ppe.map(f, segments)

    homologues = rectify_results(segment_results)

    print("Homologue collapsing runtime:", time() - start)
    return homologues

def rectify_results(segment_results):

    homologues = collections.defaultdict(set)

    to_del = set()

    for h in segment_results:
        for k in h:
            homologues[k] |= h[k]


    for key in homologues:
        to_add = set()
        for homologue in homologues[key]:
            if homologue in homologues:
                to_add |= homologues[homologue]
                to_del.add(homologue)

        homologues[key] |= to_add

    return {k: v for k, v in homologues.items() if k not in to_del}

def collapse_segment(segment, gene_dict, lengths, min_identity, min_coverage):

    to_skip = set()
    homologues = collections.defaultdict(set)

    for gene in segment:
        if gene in to_skip:
            continue

        filter_func = partial(to_search, g1=gene, lengths=lengths,
                              to_skip=to_skip, min_identity=min_identity, min_coverage=min_coverage)

        searches = list(filter(filter_func, gene_dict))

        fasta = fasta_formatter(gene_dict, *searches)

        clusters = cluster_collapse(fasta, lengths, min_identity, min_coverage)

        for cluster in clusters:
            to_skip |= clusters[cluster]
            homologues[cluster] |= clusters[cluster]

            for h in clusters[cluster]:
                homologues[cluster] |= homologues[h]
                del homologues[h]

    return homologues
