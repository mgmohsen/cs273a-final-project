import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import unittest
from utils.perm_indexer import *

class PermutationIndexerTest(unittest.TestCase):
    def test_length(self):
        # numbers directly fetched from the paper
        # Discriminative prediction of mammalian enhancers from DNA sequence
        # https://genome.cshlp.org/content/21/12/2167.long
        self.assertEqual(len(sum_reverse_complement(construct_vector(6))), 2080)
        self.assertEqual(len(sum_reverse_complement(construct_vector(7))), 8192)


if __name__ == '__main__':
    unittest.main()
