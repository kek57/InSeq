#!/usr/bin/python -w -s

## Takes the bowtie output with the parameters as
        # -m 1 -l 17 -n 1 -a --best --strata

import sys
import collections


if len(sys.argv) == 0:
    print "\n Usage: python /path/process_bowtie.py <bowtie_output_file> \n"

IN = open(sys.argv[1])

total_count = 0
count = collections.defaultdict(dict) # count[key = chromosome][data
insertion_F = collections.defaultdict(dict)
insertion_R = collections.defaultdict(dict)
sites = {}          # sites[sampleID] = number of insertion loations for that sample


for line in IN:
    if line.startswith("#"):
        data = line.strip().split('\t')
        read_ID = data[0].split(':')
        total_count += read_ID[2]               
        reads = read_ID[2]              # number of reads with that sequence in that sample
        chromosome = data[2]
        ## suppose to be some sort of count line here
        if data[1] == "+":              # Set insertion location of "T" in "TA" top-strand dinucleotide
            site = data[3] + read_ID[1] - 1
            # add read to count of forward reads for that site in that sample
            if chromosome not in insertion_F.keys():
                if site not in insertion_F[chromosome].keys():
                    insertion_F[chromosome][site] = reads
            elif site not in insertion_F[chromosome].keys():
                insertion_F[chromosome][site] = reads
            else:
                insertion_F[chromosome][site] += reads
        if data[1] == "-":              
            site = data[3] + 1
            # add read to count of reverse reads for that site in that sample
            if chromosome not in insertion_R.keys():
                if site not in insertion_R[chromosome].keys():
                    insertion_R[chromosome][site] = reads
            elif site not in insertion_F[chromosome].keys():
                insertion_R[chromosome][site] = reads
            else:
                insertion_R[chromosome][site] += reads


for chromosome in insertion_F.keys():
    sites[chromosome] = []
    outfile = str(sys.argv[1] + '_processed.txt_' + chromosome)
    OUT = open(outfile, 'w')
    for location in insertion_F[chromosome].keys():
        sites[chromosome].append(location)
        OUT.write(str(">" + chromosome + '\t' + location + '\t' + insertion_F[chromosome][location] + '\t'))
        mysum = insertion_F[chromosome][location]
        if location in insertion_R[chromosome].keys():
            OUT.write(str(insertion_R[chromosome][location] + '\t'))
            mysum = mysum + insertion_R[chromosome][location]
            del insertion_R[chromosome][location]
        else:
            OUT.write(str('0\t'))
    OUT.write(str(mysum + '\n'))

for chromosome in insertion_R.keys():
    outfile = str(sys.argv[1] + '_processed.txt_' + chromosome)
    OUT = open(outfile, 'w')
    for location in insertion_R[chromosome].keys():
        sites[chromosome].append(location)
        OUT.write(str(">" + chromosome + '\t' + location + '\t0\t' + insertion_R[chromosome][location] + '\t'))
                    
            
            
