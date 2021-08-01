import pandas as pd
import time
from Bio.Seq import Seq
from random import sample
from utils import sequence_from_guide_and_target
import RNA
import multiprocessing

N_SAMPLES = 1

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
    target = Seq(guide).reverse_complement()
    return str(target)

def get_oligos_for_guide(guide):
    return {'reverse' : str(Seq(guide + 'GGCA').reverse_complement()), 'forward' : 'CAGC' + guide}

def get_oligos_for_target(target):
    return {'reverse': str(Seq(target).reverse_complement()), 'forward': 'GGCC' + target + 'GGCC'}

"""
We want to generate 51bp guides with 'TAG' near the center 
___24p___TAG___24bp___
We'll predefine the A %; sample locations to place A, and then randomly fill in the rest with G's and C's
"""

def generate_guides(n_as=0, n_guides=10000):
    guides = []
    stop = list("TAG")
    for i in range(n_guides):
        guide = [0] * 24 + stop + [0] * 24
        sample_indices = list(range(24)) + list(range(27, 51))
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
    os_guide = generate_guides(n_as=0, n_guides=N_SAMPLES)[0]
    os_target = get_target_for_guide(os_guide)
    os_g_primers = get_oligos_for_guide(os_guide)
    os_t_primers = get_oligos_for_target(os_target)
    os = {"guide": os_guide, "target":os_target, "guide primers": os_g_primers, "target primers": os_t_primers}

    ten_guide = generate_guides(n_as=int(.1 * 50)-1, n_guides=N_SAMPLES)[0]
    ten_target = get_target_for_guide(ten_guide)
    ten_g_primers = get_oligos_for_guide(ten_guide)
    ten_t_primers = get_oligos_for_target(ten_target)
    ten = {"guide": ten_guide, "target":ten_target, "guide primers": ten_g_primers, "target primers": ten_t_primers}

    ninety_guide = generate_guides(n_as=int(.9 * 50)-1, n_guides=N_SAMPLES)[0]
    ninety_target = get_target_for_guide(ninety_guide)
    ninety_g_primers = get_oligos_for_guide(ninety_guide)
    ninety_t_primers = get_oligos_for_target(ninety_target)
    ninety = {"guide": ninety_guide, "target":ninety_target, "guide primers": ninety_g_primers, "target primers": ninety_t_primers}

    print("os", os)
    print("\nten", ten)
    print("\nninety", ninety)
    """
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
    free_energies = pd.DataFrame(list(zip(only_stop_seqs, ten_percent_seqs, ninety_percent_seqs)), columns=['two', 'ten', 'ninety'])
    free_energies.to_csv("free_energies.csv")
    print("Finished and wrote to CSV")
    """
