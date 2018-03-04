'''
Class Name: SequenceReader
Usage: sample usage in notebooks/data_pipeline.ipynb
'''

import os, sys
import numpy as np
from collections import defaultdict
# library for viewing how much processing is completed
from tqdm import tqdm
# utility functions for indexing permutation
from perm_indexer import *

class SequenceReader():
    # initialize the class with:
    #   1) target   : target (e.g. limb, hindbrain (rhombencephalon), etc)
    #   2) k        : k-mer length
    #   3) data_file: where sequence data lives at
    def __init__(self, target, k, data_file = os.path.join("../data/sequences")):
        self.target = target
        self.k = k
        self.data_file = data_file
        self.entire_features = {'True': [], 'False': []}
        self.setup_keys()

    # read the entire sequence file
    # store feature vector and metadata for each sequence based on whether
    # the sequence shows any enhancer activity to target expression or not
    def read_all(self):
        # counter for logging
        seq_count, positive_count = 0, 0
        # read the entire sequences
        f = open(self.data_file, 'rb')
        # variable for storing sequence so far (required due to file format)
        seq_so_far = ""
        # variable for storing whether current seq is positive to the target
        positive_to_target = False
        # loop over every single line in the file
        for idx, _line in enumerate(tqdm(f)):
            # strip unnecessary new line at the end for each line
            line = _line.strip()
            # if empty line, process the sequence so far
            if len(line) == 0:
                # capitalizes the sequence
                seq_so_far = seq_so_far.upper()
                # get the normalized feature vector
                feature = self.get_feature(seq_so_far)
                # store positive example to positive dict and negative to neg
                self.entire_features[str(positive_to_target)].append([feature, info])
                # reset seq_so_far and positive_to_target variables
                seq_so_far = ""
                positive_to_target = False
                continue
            # if new sequence starts
            if line.startswith('>'):
                # increment the count of sequences
                seq_count += 1
                # parse the metadata
                splits = [_split.strip() for _split in line[1:].split('|')]
                info = {}
                info['species'] = splits[0]
                chrom, interval = splits[1].split(':')
                info['chrom'] = chrom
                info['start'], info['end'] = [int(_interval) for _interval in interval.split('-')]
                info['element'] = splits[2]
                info['exists'] = splits[3]
                # if the sequence has any positive enhancer activity
                if info['exists'] == 'positive':
                    # increment the count of sequences with positive mark
                    positive_count += 1
                    info['expressions'] = {}
                    for split in splits[4:]:
                        score_idx = split.find('[')
                        expression, score = split[:score_idx], split[score_idx:]
                        info['expressions'][expression] = score
                    # if current target expression is positive, set positive_to_target to True
                    if self.target in info['expressions']:
                        positive_to_target = True
                continue
            # append current line to sequence_so_far
            seq_so_far += line
        # special case handling for the last sequence
        # (b/c there's no more empty line left, we need this handling)
        seq_so_far = seq_so_far.upper()
        feature = self.get_feature(seq_so_far)
        self.entire_features[str(positive_to_target)].append([feature, info])
        seq_so_far = ""
        positive_to_target = False
        # log the result statistics
        print "Total lines read: {}".format(idx + 1)
        print "Total sequences read: {}".format(seq_count)
        print "Number of seqs with enhancer activity: {}".format(positive_count)
        print "# of Positive sequence with {} expression : {}".format(self.target, len(self.entire_features['True']))
        print "# of Negative sequence with {} expression : {}".format(self.target, len(self.entire_features['False']))
        f.close()

    def setup_keys(self):
        _vector_dict = construct_vector(self.k)
        _vector_dict = sum_reverse_complement(_vector_dict)
        self.keys = _vector_dict.keys()

    # helper function for computing feature vector for given sequence
    def get_feature(self, seq):
        # initialize vector with dimension 4^k
        vector_dict = construct_vector(self.k)
        # for each k-mer
        for i in range(0, len(seq)-self.k+1):
            # add the count
            try:
                vector_dict[seq[i: i+self.k]] += 1
            except KeyError:
                continue
        # sum of reverse complement counts into one as mentioned in the paper
        values = []
        for key in self.keys:
            reverse_key = reverse_complement(key)
            if key == reverse_key:
                values.append(vector_dict[key])
            else:
                values.append(vector_dict[key] + vector_dict[reverse_key])
        # convert values into numpy arrays
        feature = np.array(values)
        # normalize with l2 norm -> make to unit vector
        return feature / np.linalg.norm(feature)

    # getter function for fetching the entire entire_features
    def get_entire_features(self):
        return self.keys, self.entire_features
