from functools import partial
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from numpy import percentile
from utilities import basename, contents, load_gene, match_genes
import json2csv
import os
import pandas as pd
import subprocess
import tempfile
import update_definitions
import json

def create(genome_dir, alleles_dir, json_dir, work_dir, prokka_out, min_identity, min_coverage, mist_bin, cores):
 
    call_table = os.path.join(work_dir, 'wgmlst_calls.csv')

    #annotate(genome_dir, prokka_out, cores)

    representatives = resolve_homologues(prokka_out, work_dir, min_identity, min_coverage, cores)

    create_markers(representatives, alleles_dir, work_dir)

    perform_wgmlst(mist_bin, test_path, alleles_dir, genome_dir, json_dir, cores)

    update_definitions.update_definitions(alleles_dir, json_dir, 'wgmlst')

    json2csv.convert_to_table(json_dir, call_table, work_dir)

    divide_schemes(call_table, threshold, remove_worst_genomes)

def prokka(fasta, outdir):

    name = basename(fasta)
    prokka = ('prokka',
              '--cpus', '1',
              '--outdir', os.path.join(outdir, name),
              '--locustag', name,
              '--prefix', name,
              fasta)

    subprocess.call(prokka)

def annotate(fasta_dir, output_dir, cores):
    
    with ProcessPoolExecutor(cores) as ppe:
        ppe.map(partial(prokka, outdir=output_dir),
                contents(fasta_dir, '.fasta'))

def resolve_homologues(prokka_outdir, work_dir, min_identity, min_coverage, cores):

    genes = {}

    for ffn in contents(prokka_outdir, '/*.ffn'):
        with open(ffn) as f:
            cur = dict(load_gene(g) for g in SeqIO.parse(f, 'fasta'))
            genes.update(cur)

    homologue_groups = match_genes(genes, min_identity, min_coverage, cores)

    representatives = {g: genes[g] for g in homologue_groups}

    with open(os.path.join(work_dir, 'homologues.json'), 'w') as h:
        json.dump(homologue_groups, h, sort_keys=True)

    return representatives

def create_markers(reps, alleles_dir, work_dir):

    def format_row(gene):

        return '\t'.join((gene, 'wgmlst', '1', '', '',
                          '-1', '0', gene + '.fasta', '0'))

    header = '\t'.join(('marker_name', 'test_name', 'test_type',
                        'fwd_primer', 'rev_primer', 'amplicon_size',
                        'amplicon_range', 'allelic_db_filename', 'repeat_size'))

    markers = [header]

    for gene in reps:

        markers.append(format_row(gene))

        fasta = '>1\n{}\n'.format(reps[gene])

        path = os.path.join(alleles_dir, gene + '.fasta')

        with open(path, 'w') as f:
            f.write(fasta)

    with open(os.path.join(work_dir, 'wgmlst.markers'), 'w') as t:
        t.write('\n'.join(markers))

def perform_wgmlst(mist_bin, test_path, alleles_dir,
                   genome_dir, json_dir, cores):

    def run_mist(genome):

        with tempfile.TemporaryDirectory() as d:
            cmd = (mist_bin, '-T', d, '-t', test_path, '-a', alleles_dir,
                   '-b', '-j', os.path.join(json_dir, basename(genome)), genome)

            subprocess.call(cmd)

    with ProcessPoolExecutor(cores) as p:
        p.map(run_mist, (os.path.join(genome_dir, g) for g in os.listdir(genome_dir)))

def divide_schemes(results_table, threshold, remove_worst_genomes):

    def count_bad(array):

        return sum(1 for x in array if x in (-1, 0, '?'))

    calls = pd.read_csv(results_table, header=0, index_col=0)

    genome_badness =  calls.apply(count_bad, axis=0)
    remove_thresh = percentile(genome_badness,  remove_worst_genomes)

    to_keep = [genome >= remove_thresh for genome in genome_badness]

    calls = calls[to_keep]

    gene_badness = calls.apply(count_bad, axis=1)

    cgmlst = [x / len(calls) < threshold for x in gene_badness]
    accessory = [not y for y in cgmlst]

    cgmlst_calls = calls[calls.columns[cgmlst]]
    acc_calls =  calls[calls.columns[accessory]]

    return cgmlst_calls, acc_calls

def export_schemes(outdir, work_dir, cgmlst, accessory):

    # Set up call table paths
    cg_calls = os.path.join(outdir, 'cgmlst_calls.csv')
    acc_calls = os.path.join(outdir, 'accessory_calls.csv')

    # Set up .markers file paths
    wgmlst_markers = os.path.join(work_dir, 'wgmlst.markers')
    cgmlst_markers = os.path.join(work_dir, 'cgmlst.markers')
    accessory_markers = os.path.join(work_dir, 'accessory.markers')

    # Get gene sets for each scheme
    cg_genes = set(cgmlst.columns) | {'marker_name'}
    acc_genes = set(accessory.columns) | {'marker_name'}

    # Write call tables
    cgmlst.to_csv(cg_calls)
    accessory.to_csv(acc_calls)

    # Divide wgmlst.markers into cgmlst.markers and accessory.markers
    with open(wgmlst_markers, 'r') as w:

        with open(cgmlst_markers, 'a') as c, open(accessory_markers, 'a') as a:

            for line in w:
                l = line.split()

                if l[0] in cg_genes:
                    c.write(line)

                if l[0] in acc_genes:
                    a.write(line)

