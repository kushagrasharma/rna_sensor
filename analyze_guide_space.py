import pandas as pd
import json
import time
from Bio.Seq import Seq
from Bio import SeqIO
from random import sample
from utils import trans_sequence_from_guide, iter_sample_fast
import RNA
import multiprocessing

SAMPLE_SIZE = 100

try:
    cpus = multiprocessing.cpu_count()
except NotImplementedError:
    cpus = 2   # arbitrary default

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches

def cofold_sequence(seq):
    fc = RNA.fold_compound(seq)
    _, mfe = fc.mfe_dimer()
    return mfe

def get_targets_for_sequence(seq, target_len=54):
    stop_mismatch = 'CCA'
    stop_indexes = find_all(seq, stop_mismatch)
    targets = []
    guides = []
    for stop_i in stop_indexes:
        first_possible_target_start = max(0, stop_i - (target_len - len(stop_mismatch)))
        bases_at_end = len(seq[stop_i+3:])
        last_possible_target_start = stop_i if bases_at_end >= target_len - len(stop_mismatch) else max(stop_i - (target_len - bases_at_end), 0)
        for i in range(first_possible_target_start, last_possible_target_start):
            targets.append(seq[i:i+54])
            guide = targets[-1]
            guide = list(guide)
            guide[stop_i+1-i] = 'T'
            guide = ''.join(guide)
            guide = str(Seq(guide).reverse_complement())
            guides.append(guide)
    return targets, guides

if __name__ == "__main__":
    mfe_for_gene = {}

    GFP = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"

    targets, guides = get_targets_for_sequence(GFP)
    seq_to_fold = [trans_sequence_from_guide(x) + '&' + GFP for x in guides]

    print("{} cpus for multiprocessing".format(cpus))
    # Get free energies for sequences
    pool = multiprocessing.Pool(processes=cpus)

    print("Checking GFP")

    mfe_for_gene['GFP'] = pool.map(cofold_sequence, seq_to_fold)

    print("Done with GFP")
    fasta_sequences = SeqIO.parse(open("GRCh38_latest_rna.fasta"), "fasta")

    sampled_genes = iter_sample_fast(fasta_sequences, SAMPLE_SIZE)

    for i, fasta in enumerate(sampled_genes):
        print("Checking gene {} / {}".format(i, SAMPLE_SIZE))
        gene_name, seq = fasta.id, str(fasta.seq)
        targets, guides = get_targets_for_sequence(seq)
        seq_to_fold = [trans_sequence_from_guide(x) + '&' + seq for x in guides]
        mfe_for_gene[gene_name] = pool.map(cofold_sequence, seq_to_fold)

    with open('mfe.json', 'w') as fp:
        json.dump(mfe_for_gene, fp)
