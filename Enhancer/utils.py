from Bio.Seq import reverse_complement
import itertools
import numpy as np
from tqdm import tqdm

KMER_LEN = 4
FRAME_LEN = 1500

def generate_reverse_complement_keys():
    all_keys = [''.join(item) for item in itertools.product('ATGC', repeat=KMER_LEN)]
    done = dict.fromkeys(all_keys, False)
    result = []
    for key in all_keys:
        if done[key]:
            continue
        result.append(key)
        done[key] = True
        rev = reverse_complement(key)
        done[rev] = True
    return result

def prepare_frames_features(frames):
    good_frames = [f for f in frames if f.contains_N == False]
    [f.count_kmers() for f in tqdm(good_frames)]
    return np.array([f.get_features() for f in good_frames])

def split_data(data, msk):
    X = data[msk, :136]
    y = data[msk, 136]
    return X, y
