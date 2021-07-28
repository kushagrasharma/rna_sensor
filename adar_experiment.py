import pandas as pd
import time
from Bio.Seq import Seq
from random import sample
from utils import sequence_from_guide_and_target
import RNA
import multiprocessing

try:
    cpus = multiprocessing.cpu_count()
except NotImplementedError:
    cpus = 2   # arbitrary default

def get_target_for_guide(guide):
    a_index = guide.find("TAG") + 1
    assert(a_index != -1)
    guide = list(guide)
    guide[a_index] = 'G'
    guide = ''.join(guide)
    guide = Seq(guide).reverse_complement()
    return str(guide)

def get_oligos_for_guide(guide):
    return {'reverse' : 'CAGC' + guide + 'GGCA', 'forward' : str(Seq(guide).reverse_complement())}

def get_3primeUTR_oligos_for_target(target):
    return {'reverse': str(Seq(target).reverse_complement()), 'forward': 'GGCC' + target + 'GGCC'}

"""
We want to generate 50bp guides with 'TAG' near the center 
___23bp___TAG___24bp___
We'll predefine the A %; sample locations to place A, and then randomly fill in the rest with G's and C's
"""

def generate_guides(n_as=0, n_guides=10000):
    guides = []
    stop = list("TAG")
    for i in range(n_guides):
        guide = [0] * 23 + stop + [0] * 24
        sample_indices = list(range(23)) + list(range(26, 50))
        a_indices = sample(sample_indices, n_as) 
        for i in a_indices:
            guide[i] = 'A'
        for i in range(len(guide)):
            if not guide[i]:
                guide[i] = sample(['G', 'C'], 1)[0]
        guide = ''.join(guide)
        
        guides.append(guide)
    return guides
if __name__ == "__main__":
    only_stop_guides = generate_guides(n_as=0, n_guides=10000)
    ten_percent_as_guides = generate_guides(n_as=int(.1 * 50)-1, n_guides=10000)
    ninety_percent_as_guides = generate_guides(n_as=int(.9 * 50)-1, n_guides=10000)

    only_stop_seqs = []
    ten_percent_seqs = []
    ninety_percent_seqs = []

    for i in range(len(only_stop_guides)):
        only_stop_seqs.append(sequence_from_guide_and_target(only_stop_guides[i], get_target_for_guide(only_stop_guides[i])))
        ten_percent_seqs.append(sequence_from_guide_and_target(ten_percent_as_guides[i], get_target_for_guide(ten_percent_as_guides[i])))
        ninety_percent_seqs.append(sequence_from_guide_and_target(ninety_percent_as_guides[i], get_target_for_guide(ninety_percent_as_guides[i])))

    print("{} cpus for multiprocessing".format(cpus))
    # Get free energies for sequences
    pool = multiprocessing.Pool(processes=cpus)
    only_stop_seqs = [x[1] for x in pool.map(RNA.fold, only_stop_seqs)]
    ten_percent_seqs = [x[1] for x in pool.map(RNA.fold, ten_percent_seqs)]
    ninety_percent_seqs = [x[1] for x in pool.map(RNA.fold, ninety_percent_seqs)]
