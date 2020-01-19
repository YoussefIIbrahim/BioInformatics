import numpy as np
from Enhancer.utils import prepare_frames_features, split_data, FRAME_LEN
from Enhancer.frame import Frame

chr_features_file = "chr_features.csv"

def process_chr_file(file_path, k_mers):
    with open(file_path, 'r') as f:
        content = f.read()
        seq = ''.join(content.split('\n')[1:])
        step = int(FRAME_LEN/2)
        frames = [Frame(k_mers, seq[i:i+FRAME_LEN]) for i in range(0, len(seq)-step, step)]
    return frames

def preprocess_chr(chr_frames):
    chr_features = prepare_frames_features(chr_frames)
    X_chr, _ = split_data(chr_features, np.full(len(chr_features), True))
    np.savetxt(chr_features_file, np.array(X_chr))
    return X_chr

def get_chr_data():
    return np.genfromtxt(chr_features_file)

def fill_means(y, frames):
    tmp = [p[1] for p in y]
    mean = np.mean(tmp)
    i = 0
    pred = []
    for frame in frames:
        if frame.contains_N:
            pred.append(mean)
        else:
            pred.append(tmp[i])
            i += 1
    return pred
