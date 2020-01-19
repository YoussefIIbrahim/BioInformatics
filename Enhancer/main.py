from Enhancer.utils import KMER_LEN, FRAME_LEN, generate_reverse_complement_keys, split_data
from Enhancer.preprocessor import preprocess, get_train_data
from Enhancer.chr_preprocessor import process_chr_file, preprocess_chr, fill_means, get_chr_data
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
import numpy as np

k_mers = generate_reverse_complement_keys()

preprocess(k_mers)
data = get_train_data()

msk = np.random.rand(len(data)) < 0.8
X_train, y_train = split_data(data, msk)
X_test, y_test = split_data(data, ~msk)
clf = RandomForestClassifier()
clf.fit(X_train, y_train)
y_pred = clf.predict(X_test)
print(classification_report(y_test, y_pred, labels=[0, 1]))

chr_frames = process_chr_file('chr21.fa.txt', k_mers)
preprocess_chr(chr_frames)
X_chr = get_chr_data()
y_chr = clf.predict_proba(X_chr)
chr_pred = fill_means(y_chr, chr_frames)

with open("chr21.wig", "w") as f:
    f.write("fixedStep chrom=chr21 start=0 step=750 span=1500\n")
    f.write("\n".join([str(c) for c in chr_pred]))
