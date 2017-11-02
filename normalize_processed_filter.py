#!/usr/bin/python -w -s

import sys
import glob
import collections

if len(sys.argv) < 1:
    print "\n Usage: perl $0 <input.txt>\n"

IN = open(sys.argv[1], 'r')
OUT = open(str(sys.argv[1] + '_filter_cpm.txt'), 'w')

logfile = glob.glob('*.scarf.log')
LOG = open(str(logfile[0]), 'w')

LOG.write('SampleID\tRawREadsAfterFiltered\tCoverage\tScale Factor\t\n')

data = collections.defaultdict(dict)

for line in IN:
    temp = line.strip().split('\t')
    chromosome = temp[0]
    location = temp[1]
    total_reads = temp[4]
    data[chromosome][location] = total_reads
    sites = []
    counts = []
    norm_counts = []
for chrom in data.keys():
    for loc in data[chrom].keys():
        if int(data[chrom][loc]) > 3:
            sites.append(loc)
            counts.append(data[chrom][loc])
mysum = 0
coverage = len(sites)


i = 0
while i < coverage:
    mysum += int(counts[i])
    i += 1

scale_factor = 1000000 / float(mysum)


LOG.write(str(sys.argv[1] + '\t' + str(mysum) + '\t' + str(coverage) + '\t' + str(scale_factor) + '\n'))

for chrom in data.keys():
    for loc in data[chrom].keys():
        if int(data[chrom][loc]) > 3:
            norm_data = int(data[chrom][loc]) * scale_factor
            OUT.write(str(str(chrom) + '\t' + str(loc) + '\t' + str(norm_data) + '\n'))


LOG.close()
OUT.close()


    
            

