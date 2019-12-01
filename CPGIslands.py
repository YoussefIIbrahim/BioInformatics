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

    prevChrom = data.iloc['chromosome']
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
with open("data.bed")as f:
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
fig, ax = plt.subplots()
plt.bar(['sea','cpg island','shore','shelve'],hits)
plt.title('DNA Methylation hits per region')
plt.ylabel('Percentage of hits')
plt.xlabel('Region')
plt.show()
np.sum(pd.Series(regions).value_counts())




#
# #%% find islands, shores, shelves and seas: function definitions
#
# def add_to_interval(interval, start, stop, dilation, chrom_borders):
#     interval.add(spans.intrange(start - dilation, stop + dilation).intersection(chrom_borders))
#
# def intervals_for_chrom(islands_chrom_df, chrom_len):
#     shore_limit = 2000
#     shelf_limit = shore_limit + 2000
#     chrom_interval = spans.intrange(0, chrom_len)
#     intervals = {SEA_KEY: spans.intrangeset([chrom_interval]), SHELF_KEY: spans.intrangeset([]),
#                  SHORE_KEY: spans.intrangeset([]), ISLAND_KEY: spans.intrangeset([])}
#     for _, island in islands_chrom_df.iterrows():
#         add_to_interval(intervals[ISLAND_KEY], island[START_KEY], island[STOP_KEY], 0, chrom_interval)
#         add_to_interval(intervals[SHORE_KEY], island[START_KEY], island[STOP_KEY], shore_limit, chrom_interval)
#         add_to_interval(intervals[SHELF_KEY], island[START_KEY], island[STOP_KEY], shelf_limit, chrom_interval)
#     intervals[SEA_KEY] = intervals[SEA_KEY].difference(intervals[SHELF_KEY])
#     intervals[SHELF_KEY] = intervals[SHELF_KEY].difference(intervals[SHORE_KEY])
#     intervals[SHORE_KEY] = intervals[SHORE_KEY].difference(intervals[ISLAND_KEY])
#     return intervals
#
#
# def df_from_interval(interval, chrom):
#     fragments_list =[]
#     for subinterval in interval:
#         fragments_list.append((chrom, subinterval.lower, subinterval.upper))
#     df = pd.DataFrame(fragments_list)
#     return df
#
# #%% find islands, shores, shelves and seas
#
#
# region_chrom_dfs = {ISLAND_KEY: [], SHORE_KEY: [], SHELF_KEY: [], SEA_KEY: []}
# chrom_interval_dict = {}
# for chrom in autosomal_chrom:
#     islands_chrom_df = islands_df[islands_df[CHROM_KEY] == chrom]
#     intervals = intervals_for_chrom(islands_chrom_df, lengths[chrom])
#     for key in intervals:
#         region_chrom_dfs[key].append(df_from_interval(intervals[key], chrom))
#     chrom_interval_dict[chrom] = intervals
#
# region_dfs = {}
# for region in region_chrom_dfs:
#     region_dfs[region] = pd.concat(region_chrom_dfs[region])
#
# #%% save results to files
# for region in region_dfs:
#     region_dfs[region].to_csv(f'{region}.bed', sep='\t', header=False, index=False)
#
# #%% calculate methylations positions on chromosomes
# methyalation_df[POS_KEY] = (methyalation_df[STOP_KEY] - methyalation_df[START_KEY]) / 2 + methyalation_df[START_KEY]
# methyalation_df[POS_KEY] = methyalation_df[POS_KEY].apply(lambda x: int(x))
#
# #%% calculate methylations locations in islands/shores/shelves/seas
# def detemine_methylation_location(chrom, pos, chrom_interval_dict):
#     intervals = chrom_interval_dict[chrom]
#     for region in intervals:
#         if intervals[region].contains(pos):
#             return region
#     raise LookupError(f'There is no position {str(pos)} in chromosome {chrom}.')
#
#
# methyalation_df[LOCATION_KEY] = methyalation_df.apply(
#     lambda row: detemine_methylation_location(row[CHROM_KEY], row[POS_KEY], chrom_interval_dict), axis=1)
#
# #%% collect statistics about methylations locations
#
# def count_methylations_in_area(methyalation_df, area):
#     return sum(methyalation_df[LOCATION_KEY] == area)
#
#
# areas = (ISLAND_KEY, SHORE_KEY, SHELF_KEY, SEA_KEY)
# count = []
# for area in areas:
#     count.append(count_methylations_in_area(methyalation_df, area))
#
# #%% plot results
#
# def autolabel(rects):
#     """Attach a text label above each bar in *rects*, displaying its height."""
#     for rect in rects:
#         height = rect.get_height()
#         ax.annotate('{}'.format(height),
#                     xy=(rect.get_x() + rect.get_width() / 2, height),
#                     xytext=(0, 3),  # 3 points vertical offset
#                     textcoords="offset points",
#                     ha='center', va='bottom')
#
#
# x = np.arange(len(areas))  # the label locations
# width = 0.6  # the width of the bars
#
# fig, ax = plt.subplots()
# rects = ax.bar(x, count, width)
#
# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Number of methylations')
# ax.set_title('Number of methylations in each area')
# ax.set_xticks(x)
# ax.set_xticklabels(areas)
# # ax.legend()
#
# autolabel(rects)
#
# fig.tight_layout()
#
# plt.show()
#
