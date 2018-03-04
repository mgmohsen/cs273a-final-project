import itertools
import math

LETTER_DICT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
NUM_LETTERS = 4

# function that creates a dictionary for 4^k k-mer vector with 0 counts
def construct_vector(k):
    vector_dict = {}
    for key in [''.join(p) for p in itertools.product(LETTER_DICT.keys(), repeat=k)]:
        vector_dict[key] = 0
    return vector_dict

def sum_reverse_complement(vector_dict):
    # helper function for retrieving reverse complement of a given DNA
    def reverse_complement(dna):
        return ''.join([COMPLEMENT[base] for base in dna[::-1]])
    # for all 4^k keys
    keys = vector_dict.keys()
    for key in keys:
        # get the reverse complement key
        reverse_key = reverse_complement(key)
        # if itself is the reverse complement, skip this process
        if key == reverse_key:
            continue
        # if both original and reverse complement is in the dictionary
        if key in vector_dict and reverse_key in vector_dict:
            # sum two up and revmoe the latter from vector_dict
            vector_dict[key] += vector_dict[reverse_key]
            del vector_dict[reverse_key]
    return vector_dict
