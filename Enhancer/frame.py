from Bio.Seq import reverse_complement
import itertools
import numpy as np
from Enhancer.utils import FRAME_LEN, KMER_LEN

class Frame:
    def __init__(self, k_mers, seq, label=None):
        self.sequence = seq.upper()
        self.counter = dict.fromkeys(k_mers, 0)
        self.label = label
        self.contains_N = self.sequence.find('N') != -1

    def count_kmers(self):
        for i in range(len(self.sequence) - KMER_LEN + 1):
            k_mer = self.sequence[i:i + KMER_LEN]
            if k_mer in self.counter:
                self.counter[k_mer] += 1
            else:
                rev = reverse_complement(k_mer)
                self.counter[rev] += 1

    def get_features(self):
        return np.array([self.counter[key]/FRAME_LEN for key in sorted(self.counter.keys())] + [self.label])
