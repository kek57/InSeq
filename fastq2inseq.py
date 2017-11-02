#!usr/bin/python -w -s

import sys
import glob
import gzip
import os



def fastq2inseq(allfastqgz,outfile):

    OUT = open(outfile, 'w')
    
    for file in allfastqgz:                         
        IN = gzip.open(file, 'r')

        for line in IN:
            if line.startswith('@'):
                newline = '>' + line.strip()
                seq = IN.next()
                newline = newline + ':' + seq.strip()
            if line.startswith('+'):
                nextline = IN.next()
                newline = newline + ':' + nextline.strip() + '\n'
                OUT.write(newline)

        IN.close()
    OUT.close()


###################################################################
fastqpath = str(sys.argv[2] + '*.fastq*')
allfastqgz = glob.glob(fastqpath)
outfile = str(sys.argv[2] + sys.argv[1] + 'pipeline_input.txt')

fastq2inseq(allfastqgz, outfile)

    
