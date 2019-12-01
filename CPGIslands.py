import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt


def parse_chromosome(pandas_series):
    return pandas_series.apply(lambda x: int(x[3:]))

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
cpg_islands['chromosome'] = parse_chromosome(cpg_islands['chromosome'])
cpg_islands = cpg_islands.sort_values(['chromosome', 'start'])

chromosomes_lengths = pd.read_csv('hg19.chrom.sizes.txt', sep = '\t', header = None, names = ['chromosome','length'])
chromosomes_lengths = chromosomes_lengths[chromosomes_lengths['chromosome'].isin(chromosomes)]
chromosomes_lengths['chromosome'] = parse_chromosome(chromosomes_lengths['chromosome'])
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

#%% find islands, shores, shelves and seas: function definitions

def add_to_interval(interval, start, stop, dilation, chrom_borders):
    interval.add(spans.intrange(start - dilation, stop + dilation).intersection(chrom_borders))

def intervals_for_chrom(islands_chrom_df, chrom_len):
    shore_limit = 2000
    shelf_limit = shore_limit + 2000
    chrom_interval = spans.intrange(0, chrom_len)
    intervals = {SEA_KEY: spans.intrangeset([chrom_interval]), SHELF_KEY: spans.intrangeset([]),
                 SHORE_KEY: spans.intrangeset([]), ISLAND_KEY: spans.intrangeset([])}
    for _, island in islands_chrom_df.iterrows():
        add_to_interval(intervals[ISLAND_KEY], island[START_KEY], island[STOP_KEY], 0, chrom_interval)
        add_to_interval(intervals[SHORE_KEY], island[START_KEY], island[STOP_KEY], shore_limit, chrom_interval)
        add_to_interval(intervals[SHELF_KEY], island[START_KEY], island[STOP_KEY], shelf_limit, chrom_interval)
    intervals[SEA_KEY] = intervals[SEA_KEY].difference(intervals[SHELF_KEY])
    intervals[SHELF_KEY] = intervals[SHELF_KEY].difference(intervals[SHORE_KEY])
    intervals[SHORE_KEY] = intervals[SHORE_KEY].difference(intervals[ISLAND_KEY])
    return intervals


def df_from_interval(interval, chrom):
    fragments_list =[]
    for subinterval in interval:
        fragments_list.append((chrom, subinterval.lower, subinterval.upper))
    df = pd.DataFrame(fragments_list)
    return df

#%% find islands, shores, shelves and seas


region_chrom_dfs = {ISLAND_KEY: [], SHORE_KEY: [], SHELF_KEY: [], SEA_KEY: []}
chrom_interval_dict = {}
for chrom in autosomal_chrom:
    islands_chrom_df = islands_df[islands_df[CHROM_KEY] == chrom]
    intervals = intervals_for_chrom(islands_chrom_df, lengths[chrom])
    for key in intervals:
        region_chrom_dfs[key].append(df_from_interval(intervals[key], chrom))
    chrom_interval_dict[chrom] = intervals

region_dfs = {}
for region in region_chrom_dfs:
    region_dfs[region] = pd.concat(region_chrom_dfs[region])

#%% save results to files
for region in region_dfs:
    region_dfs[region].to_csv(f'{region}.bed', sep='\t', header=False, index=False)

#%% calculate methylations positions on chromosomes
methyalation_df[POS_KEY] = (methyalation_df[STOP_KEY] - methyalation_df[START_KEY]) / 2 + methyalation_df[START_KEY]
methyalation_df[POS_KEY] = methyalation_df[POS_KEY].apply(lambda x: int(x))

#%% calculate methylations locations in islands/shores/shelves/seas
def detemine_methylation_location(chrom, pos, chrom_interval_dict):
    intervals = chrom_interval_dict[chrom]
    for region in intervals:
        if intervals[region].contains(pos):
            return region
    raise LookupError(f'There is no position {str(pos)} in chromosome {chrom}.')


methyalation_df[LOCATION_KEY] = methyalation_df.apply(
    lambda row: detemine_methylation_location(row[CHROM_KEY], row[POS_KEY], chrom_interval_dict), axis=1)

#%% collect statistics about methylations locations

def count_methylations_in_area(methyalation_df, area):
    return sum(methyalation_df[LOCATION_KEY] == area)


areas = (ISLAND_KEY, SHORE_KEY, SHELF_KEY, SEA_KEY)
count = []
for area in areas:
    count.append(count_methylations_in_area(methyalation_df, area))

#%% plot results

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


x = np.arange(len(areas))  # the label locations
width = 0.6  # the width of the bars

fig, ax = plt.subplots()
rects = ax.bar(x, count, width)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of methylations')
ax.set_title('Number of methylations in each area')
ax.set_xticks(x)
ax.set_xticklabels(areas)
# ax.legend()

autolabel(rects)

fig.tight_layout()

plt.show()

