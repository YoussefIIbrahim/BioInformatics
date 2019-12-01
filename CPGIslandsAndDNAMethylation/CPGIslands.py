import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

#%% definitions
LOCATION_KEY = 'location'
POS_KEY = 'pos'
CHROM_KEY = 'chr'
STOP_KEY = 'stop'
START_KEY = 'start'
SEA_KEY = 'sea'
SHELF_KEY = 'shelf'
SHORE_KEY = 'shore'
ISLAND_KEY = 'island'


def parseChromosome(pandas_series):
    return pandas_series.apply(lambda x: int(x[3:]))


def saveBedFormat(filename, data):
    with open(filename+'.bed', 'w') as f:
        for row in data:
            f.write("%s\t%s\t%s\n" % row)


def constructCoordinatesBed(data, window=2000):
    CPGIslands = []
    shores = []
    shelves = []
    seas = []

    prevChrom = data.iloc[0]['chromosome']
    prevStart = data.iloc[0]['start']
    prevStop = data.iloc[0]['stop']
    prevChromLen = data.iloc[0]['length']

    CPGIslands.append((prevChrom, prevStart, prevStop))
    if prevStart > (2*window):
        shores.append((prevChrom, prevStart - window, prevStart))
        shelves.append((prevChrom, prevStart - 2*window, prevStart - window))
        seas.append((prevChrom, 0, prevStart - 2*window))
    elif prevStart > window:
        shores.append((prevChrom, prevStart - window, prevStart))
        shelves.append((prevChrom, 0, prevStart - window))
    elif prevStart < window:
        shores.append((prevChrom, 0, prevStart))
    else:
        sys.exit("Error: Window value is not compatible with init file structure")

    for i in tqdm(range(1, data.shape[0]), desc = 'CPG Islands'):
        currStart = data.iloc[i]['start']
        currStop = data.iloc[i]['stop']
        currChrom = data.iloc[i]['chromosome']
        currChromLen = data.iloc[i]['length']

        if prevChrom == currChrom:
            if currStart - prevStop > (4 * window):
                shores.append((currChrom, prevStop, prevStop + window))
                shelves.append((currChrom, prevStop + window, prevStop + 2 * window))
                seas.append((currChrom, prevStop + 2 * window, currStart - 2 * window))
                shelves.append((currChrom, currStart - 2 * window, currStart - window))
                shores.append((currChrom, currStart - window, currStart))
            elif currStart - prevStop > (2 * window):
                shores.append((currChrom, prevStop, prevStop + window))
                shelves.append((currChrom, prevStop + window, currStart - window))
                shores.append((currChrom, currStart - window, currStart))
            elif currStart - prevStop > window or currStart - prevStop <= window:
                shores.append((currChrom, prevStop, currStart))
            else:
                sys.exit("Something's wrong")
        else:
            if prevChromLen - prevStop > (2 * window):
                shores.append((prevChrom, prevStop, prevStop + window))
                shelves.append((prevChrom, prevStop + window, prevStop + 2 * window))
                seas.append((prevChrom, prevStop + 2 * window, prevChromLen))
            elif prevChromLen - prevStop > (window):
                shores.append((prevChrom, prevStop, prevStop + window))
                shelves.append((prevChrom, prevStop + window, prevChromLen))
            elif prevChromLen - prevStop <= (window):
                shores.append((prevChrom, prevStop, prevChromLen))
            else:
                sys.exit("Something's wrong")

            # CURRENT side
            if currStart > (2 * window):
                seas.append((currChrom, 0, currStart - 2 * window))
                shelves.append((currChrom, currStart - 2 * window, currStart - window))
                shores.append((currChrom, currStart - window, currStart))
            elif currStart > (window):
                shelves.append((currChrom, 0, currStart - window))
                shores.append((currChrom, currStart - window, currStart))
            elif currStart <= (window):
                shores.append((currChrom, 0, currStart))
            else:
                sys.exit("Something's wrong")

        prevStart = currStart
        prevStop = currStop
        prevChrom = currChrom
        prevChromLen = currChromLen
        CPGIslands.append((currChrom, currStart, currStop))

    if prevChromLen - prevStop > (2 * window):
        shores.append((prevChrom, prevStop, prevStop + window))
        shelves.append((prevChrom, prevStop + window, prevStop + 2 * window))
        seas.append((prevChrom, prevStop + 2 * window, prevChromLen))
    elif prevChromLen - prevStop > window:
        shores.append((prevChrom, prevStop, prevStop + window))
        shelves.append((prevChrom, prevStop + window, prevChromLen))
    elif prevChromLen - prevStop <= window:
        shores.append((prevChrom, prevStop, prevChromLen))
    else:
        sys.exit("Something's wrong")

    return CPGIslands, shelves, shores, seas


def calculate_hits(CPG_Islands, shores, shelves, seas, DNA_Methylation):
    """
    All data must be in tuples
    (chromosome, start, stop)

    """
    # Prepare data
    cpg_islands_df = pd.DataFrame(CPG_Islands, columns=['chromosome', 'start', 'stop'])
    cpg_islands_df['region'] = 'cpg_island'
    shores_df = pd.DataFrame(shores, columns=['chromosome', 'start', 'stop'])
    shores_df['region'] = 'shore'
    shelves_df = pd.DataFrame(shelves, columns=['chromosome', 'start', 'stop'])
    shelves_df['region'] = 'shelve'
    seas_df = pd.DataFrame(seas, columns=['chromosome', 'start', 'stop'])
    seas_df['region'] = 'sea'
    all_data = cpg_islands_df.append(shores_df, ignore_index=True).append(shelves_df, ignore_index=True).append(seas_df,
                                                                                                                ignore_index=True)
    all_data = all_data.sort_values(['chromosome', 'start'])

    regions = []
    for row in tqdm(DNA_Methylation, desc='DNA Methylation'):
        middle = row[1] + (row[2] - row[1]) / 2
        all_data_rows = all_data[
            (all_data['chromosome'] == row[0]) & (all_data['start'] <= middle) & (all_data['stop'] >= middle)]

        if all_data_rows.shape[0] > 1:
            all_data_rows = all_data_rows.iloc[[-1]]
        if all_data_rows.shape[0] == 0:
            print('no match')
            print(row)

        regions.append(all_data_rows['region'].values)
    return regions




#%% read files
# islands_df = pd.read_csv('cpgIslandExt.txt', sep='\t', usecols=[1, 2, 3], header=None)
# islands_df.columns=[CHROM_KEY, START_KEY, STOP_KEY]
# lengths_df = pd.read_csv('hg19.chrom.sizes.txt', sep='\t', index_col=0, header=None, squeeze=True)
# lengths = lengths_df.to_dict()
# methyalation_df = pd.read_csv('data.bed', sep='\t', usecols=[0, 1, 2], header=None)
# methyalation_df.columns=[CHROM_KEY, START_KEY, STOP_KEY]

cpg_islands = pd.read_csv('cpgIslandExt.txt', sep = '\t', header = None).iloc[:,1:4]
cpg_islands = cpg_islands.rename(columns = {1:'chromosome',2:'start',3:'stop'})
chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14',
               'chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
cpg_islands = cpg_islands[cpg_islands['chromosome'].isin(chromosomes)]
cpg_islands['chromosome'] = parseChromosome(cpg_islands['chromosome'])
cpg_islands = cpg_islands.sort_values(['chromosome', 'start'])

chromosomes_lengths = pd.read_csv('hg19.chrom.sizes.txt', sep = '\t', header = None, names = ['chromosome','length'])
chromosomes_lengths = chromosomes_lengths[chromosomes_lengths['chromosome'].isin(chromosomes)]
chromosomes_lengths['chromosome'] = parseChromosome(chromosomes_lengths['chromosome'])
data = cpg_islands.merge(chromosomes_lengths)

dna_methylation = []
with open("hg38_CpG-Island.bed")as f:
    for line in f:
        row = line.strip().split()
        if row[0] in chromosomes:
            chromosome = int(row[0][3:])
            start = int(row[1])
            stop = int(row[2])
            dna_methylation.append((chromosome, start, stop))

cpg_islands, shelves, shores, seas = constructCoordinatesBed(data, 2000)
saveBedFormat('cpg_islands', cpg_islands)
saveBedFormat('shelves', shelves)
saveBedFormat('shores', shores)
saveBedFormat('seas', seas)

regions = calculate_hits(cpg_islands, shores, shelves, seas, dna_methylation)
hits = pd.Series(regions).value_counts(normalize=True)
# fig, ax = plt.subplots()
# plt.bar(['sea','cpg island','shore','shelve'],hits)
# plt.title('DNA Methylation hits per region')
# plt.ylabel('Percentage of hits')
# plt.xlabel('Region')
# plt.show()
np.sum(pd.Series(regions).value_counts())
