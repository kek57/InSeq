#!/usr/bin/python -w -s

import os
import sys
import getopt
import os.path
import collections
import glob
import re
import time

start_time = time.time()

if len(sys.argv) < 6 or len(sys.argv) > 14:
    print "Usage: perl $0 -i <the raw reads file>  -m <Barcodes mapping file>  -s <indexed genome name>  -d <length_disrupt_percent (max=1)> [-operon -c <operon_probability_cutoff (max=1)>] [-arrayed]\n"
    print "Required argument: "
    print "-i gives the input raw reads file, -m gives the mapping file with the barcode sequence and name for each sample in the format as <barcode>\t<Sample name>, -s gives the name of the indexed genome that the reads should map to \n"
    print "Optional argument: \n"
    print "-d gives the region of the gene in which insertions are expected to disrupt gene function. The default is 1, which means when insertion falls anywhere in the gene (100%), the gene function will be disrupted. Setting the -d argument to 0.9, for example, would exclude insertions in the distal 10% of the gene when calculating the total number of reads/insertions for that gene. -operon (no argument) specifies that putative downstream (polar) insertions should be calculated based on a user-provided operon probability file. -c is the cutoff for operon probability, default is 1, which means only when the probability of two genes being in an operon is equal to 1 (100%), the polar effect will be considered for the downstream gene. Setting the -c argument to 0.8, for example, will calculate a polar effect for the downstream gene if the probability of the genes being in an operon is at least 0.8 (80%). -arrayed (no argument) is the option for the arrayed library. Refer to README for more information\n"

datasource = ""                   #### The raw reads file which will be analyzed
mismatch1 = 1                     #### Mismatches allowed in finding the transposon sequence
mismatch2 = 2                     #### Mismatches allowed when mapping the sequences to reference using bowtie
mapping = ""                      #### Mapping file name
poolname = "INSEQ_experiment"     #### The prefix of all the analyzed files
index = ""                        #### The name of the index fold in the indexes
operon = 0                        #### The option if analyze the sequence using operon information
cutoff = 100                                                #### The cutoff for operon probability, default is 100
disruption = 100                                            #### The disruptable percentage of genes. Only when transposon was inserted into this proximal region of gene, the gene is considered interrupted since the insertion in the distal region of a gene may not affect the function of gene at all. The default is 100.
arrayed = 0                                                 #### The option if the data is from an arrayed library
scriptpath = os.path.abspath(os.path.dirname(sys.argv[0]))  #### Get the path for this analysis package
indexpath = os.path.abspath(os.path.dirname(sys.argv[6]))   #### Get the path of the indexes         
outpath = os.getcwd().rstrip()                              #### Get the current working directory for output

options, args = getopt.getopt(sys.argv[1:], "i:s:m:o:d:c:a", ["input=","index=","mappingfile=","operon","disruption=","cutoff=","arrayed"])

for opt in options:
    if opt[0] in ('-i', '--input'):
        datasource = opt[1]
    elif opt[0] in ('-s', '--species'):
        index = opt[1]
    elif opt[0] in ('-m', '--mappingfile'):
        mapping = opt[1]
    elif opt[0] in ('-o', '--operon'):
        operon = 1
    elif opt[0] in ('-d', '--disruption'):
        disruption  = opt[1]
    elif opt[0] in ('-c', '--cutoff'):
        cutoff = opt[1]
    elif opt[0] in ('-a', '--arrayed'):
        arrayed = 1
    else:
        print "Too many arguments"

bowtie_path = ""
try:
    config = open(str(scriptpath + '/config.txt'),'r')
except IOError:
    print "Error, please check configuration file as config.txt in the package"
for line in config:
    if re.match('bowtie_dir="/', line) and re.match(r'.*/.*/.*', line):
        line_parsed = line.split('bowtie_dir=')
        bowtie_path = line_parsed[1][1:-2]
        #print bowtie_path
        
    elif re.match('bowtie_dir="/', line) and not re.match(r'.*/.*/.*', line):
        print "please specify where the bowtie directory is\n"

outputdir = os.getcwd().rstrip()



myinput = str(poolname + ".scarf")
if os.path.isfile(myinput) == True:
    os.system(str('rm ' + myinput + '\n'))
os.system(str('ln -s ' + datasource + ' ' + myinput + '\n'))    # makes symbolic link in cwd to raw reads

if os.path.isdir('bcsortedseq') == True:
    os.system(str('rm ./bcsortedseq/*\n'))
    os.system(str('rmdir ./bcsortedseq\n'))
os.mkdir('bcsortedseq')

if os.path.isdir('results') == True:
    os.system(str('rm ./results/*\n'))
    os.system(str('rmdir ./results\n'))
os.mkdir('results')

elapsed = str(round((time.time() - start_time)/60))

print str("@ " + elapsed + " minutes -- Inputs correct!")

### Step 1 ##################################################
##### Use the mapping file to assign each read to a barcode and store the output file as inputfile_assigned.txt",
    ##### and store the statistics in the log file

barcode_assigned = str(myinput + '_assigned.txt')
logfile = str(myinput + '.log')

if os.path.isfile(mapping):
    codes = open(mapping, 'r') # file containing barcode with associated sample name
else:
    print " Error: can't open the mapping file, check the README for more details about mapping file\n"

codes_dictionary = {}       # codes_dictionary[code] = sampleID
codes_number = {}           # codes_number[code] = sample number (within a given sequencing lane)
sample_count = 1
has_barcode = {}            # has_barcode[code] = number of reads with barcode code
total = 0
total_mapped = 0
code_length = 0
                                  
for line in codes:
    codes_array = line.strip().split()
    code_length = len(codes_array[0])
    codes_dictionary[codes_array[0]] = codes_array[1]
    codes_number[codes_array[0]] = sample_count
    sample_count += 1
                                  
codes.close()

## To a single otuput file, writes lines
    # >SampleName:Barcode
    # the rest of the sequence associated with that barcode
OUT = open(barcode_assigned, 'w')
INPUT = open(myinput, 'r')

for line in INPUT:
    temp_array = line.strip().split(':')
    total += 1
    seq_barcode = temp_array[10][0:code_length]
    if seq_barcode in codes_dictionary.keys():          # If the first bases of the read are one of the barcodes
        
        total_mapped += 1
        if seq_barcode not in has_barcode.keys():
            has_barcode[seq_barcode] = 1
        else:
            has_barcode[seq_barcode] += 1
        seq = temp_array[10][code_length:]
        
        OUT.write(str('>' + codes_dictionary[seq_barcode] + ":" + seq_barcode + "\n" + seq + "\n"))
OUT.close()

elapsed = str(round((time.time() - start_time)/60))

print str("@ " + elapsed + "minutes -- Done assigning")

## Creates the appropriate files for output of trimmed sequences
filehandledictionary = {}
for key in codes_dictionary:
    if key in has_barcode.keys():
        filehandle = str(outpath + '/bcsortedseq/' + myinput + '_' + codes_dictionary[key] + '.fas')
        filehandledictionary[key] = filehandle

elapsed = str(round((time.time() - start_time)/60))

print str("@ " + elapsed +  "minutes -- File handles made")
                                  
## Writes some basic statists to LOG file
LOG = open(logfile, 'w')
LOG.write(str("Number of total reads in " + myinput + " is " + str(total) + "\n"))
LOG.write(str("Total mapped " + str(total_mapped) + "\t" + str(100 * total_mapped / total) + " \n"))
LOG.write(str("Barcode\tSample\tReads\tPercent\n"))

for key in codes_dictionary:
    barcode = key
    sample = codes_dictionary[key]
    reads = has_barcode[key]
    percent = 100 * reads / total
    LOG.write(str(barcode + '\t' + sample + '\t' + str(reads) + '\t' + str(percent) + '\n'))

elapsed = str(round((time.time() - start_time)/60))

print str("@ " + elapsed + " minutes -- Log written, beginning trimming!")

### Step 2 ##################################################
##### Trimmed reads to remove transposon, append 16bp reads with a 5'N    

## Opens previously made assigned barcode file to parse through lines and trim sequences
    
IN = open(barcode_assigned, 'r')              

tn_array = 'ACAGGTTG'
total = 0
count = 0
exact_match_count = 0
mismatch1_count = 0
mismatch2_count = 0
tn_match = {} # tn_match[sampleID] = number of reads that have TN sequence at allowable mismatches
header = []
bc = ''
number = 0
length = 0
number_dictionary = collections.defaultdict(dict) # number_dictionary[key = barcode][innerkey = trimmed_sequence] = number
length_dictionary = collections.defaultdict(dict) # number_dictionary[key = barcode][innerkey = trimmed_sequence] = length
tot = 0
for line in IN:
    line = line.strip()
    tot += 1
    print tot
# if the line is a header line do this:
    if line.startswith('>'):
        parsed_line = line.split(':')
        bc = parsed_line[1]
	print bc
        if parsed_line[0] not in tn_match.keys():
            tn_match[parsed_line[0]] = 0
# if the line is a sequence line do this:
    else:
        seq = line.strip()
        pos1_match = 0
        pos2_match = 0
    # does read contain transposon at bp 20 or 21
        tn_pos1 = seq[16:24]
        tn_pos2 = seq[17:25]

        if tn_pos1 == 'ACAGGTTG':
            pos1_match = 1
            exact_match_count += 1
            count +=1
            
        if tn_pos2 == 'ACAGGTTG':
            pos2_match = 1
            exact_match_count += 1
            count += 1
            
        if mismatch1 == 1 and tn_pos1 != 'ACAGGTTG' and tn_pos2 != 'ACAGGTTG':
            pos1_score = 0
            pos2_score = 0
            for i in range(1,len(tn_pos1)):
                if tn_pos1[i] == tn_array[i]:
                    pos1_score += 1
            for i in range(1,len(tn_pos2)):
                if tn_pos2[i] == tn_array[i]:
                    pos2_score += 1
            if pos1_score == 7:
                pos1_match = 1
                mismatch1_count += 1
                count += 1
            if pos2_score == 7:
                pos2_match = 1
                mismatch1_count += 1
                count += 1
                
        if mismatch1 == 2 and tn_pos1 != 'ACAGGTTG' and tn_pos2 != 'ACAGGTTG':
            pos1_score = 0
            pos2_score = 0
            for i in range(1,len(tn_pos1)):
                if tn_pos1[i] == tn_array[i]:
                    pos1_score += 1
            for i in range(1,len(tn_pos2)):
                if tn_pos2[i] == tn_array[i]:
                    pos2_score += 1
            if pos1_score == 7:
                pos1_match = 1
                mismatch1_count += 1
                count += 1
            if pos2_score == 7:
                pos2_match = 1
                mismatch1_count += 1
                count += 1
            if pos1_score == 6:
                pos1_match = 1
                mismatch2_count += 1
                count += 1
            if pos2_score == 6:
                pos2_match = 1
                mismatch2_count += 1
                count += 2
    # if match found

        if pos1_match == 1:
            trimmed_seq = seq[0:16]    # trim transposon sequence
            tn_match[parsed_line[0]] += 1

            if bc not in number_dictionary.keys():
                if trimmed_seq not in number_dictionary[bc].keys():
                    number_dictionary[bc][trimmed_seq] = 1
            elif trimmed_seq not in number_dictionary[bc].keys():
                number_dictionary[bc][trimmed_seq] = 1
            else:
                number_dictionary[bc][trimmed_seq] += 1
          
            length_dictionary[bc][trimmed_seq] = 16
            #print str('pos1 match found! ' + str(number_dictionary[bc][trimmed_seq]))

        if pos2_match == 1:
            trimmed_seq = seq[0:17]
            tn_match[parsed_line[0]] += 1

            if bc not in number_dictionary.keys():
                if trimmed_seq not in number_dictionary[bc].keys():
                    number_dictionary[bc][trimmed_seq] = 1
            elif trimmed_seq not in number_dictionary[bc].keys():
                number_dictionary[bc][trimmed_seq] = 1
            else:
                number_dictionary[bc][trimmed_seq] += 1
          
            length_dictionary[bc][trimmed_seq] = 17
            #print str('pos2 match found! ' + str(number_dictionary[bc][trimmed_seq]))

IN.close()

elapsed = str(round((time.time() - start_time)/60))

print str("@ " + elapsed +  " minutes -- Done trimming!")

# Export trimmed sequences

for sample in number_dictionary.keys():
    fh = filehandledictionary[sample]
    OUT = open(fh, 'w')
    for read in number_dictionary[sample].keys():
        header = str('>' + sample + ':' + str(length_dictionary[sample][read]) + ':' + str(number_dictionary[sample][read]) + '\n')
        OUT.write(header)
        OUT.write(str(read) + '\n')
    OUT.close()


LOG.write('In the trimming process, here is the statistics for each sample\n')
LOG.write('Sample\tTrimmed\tPercentage\n')
for key in codes_number.keys():
    if has_barcode[key]:
        trimkey = str('>' + codes_dictionary[key])
        percentage = str(round(tn_match[trimkey] / float(has_barcode[key]) * 100,2))
        LOG.write(str(key) + '\t' + str(tn_match[trimkey]) + '\t' + str(percentage) + '\n')
LOG.close()

os.chdir(str(index)) ### The reference genome/indexes need to be in a folder called indexes in the cwd
ptt = glob.glob('*.ptt')

os.chdir(str(outputdir))

clean_up = 'clean_up.sh'
CL = open(clean_up, 'wr')
CL.write('rm barcode_assigned \n rm mappingjobs_*\n rm -r bcsortedseq \n')

elapsed = str(round((time.time() - start_time)/60))

print str("@ " + elapsed +  " minutes -- Done exporting!")

# Create job files for each sample

elapsed = str(round((time.time() - start_time)/60))

print str("@ " + elapsed +  " minutes -- Creating jobs!")

for key in codes_dictionary.keys():
    if key in has_barcode.keys():
        path_results = str(outpath + '/results/' + myinput +'_' + codes_dictionary[key] + '.bowtiemap')
        path_indexes = str(index)
        ## Create job file
        TASKSFORBC = open(outputdir + '/mappingjobs_' + myinput + '_' + key + '.job', 'w')
        ## Print intro to job file
        TASKSFORBC.write(str("#!/bin/bash \n #SBATCH --partition=general \n #SBATCH --job-name=my_job \n #SBATCH --ntasks=1 --nodes=1 \n #SBATCH --mem-per-cpu=6000 \n #SBATCH --time=96:00:00 \n #SBATCH --mail-type=ALL \n #SBATCH --mail-user=kaitlyn.kortright@yale.edu\n"))
        ## Print bowtie command to job file
        TASKSFORBC.write(str(bowtie_path + '/bowtie -m 1 --best --strata -a --fullref -n mismatch2 -l 17 ' + path_indexes + ' -f ' + outpath + '/bcsortedseq/' + myinput + '_' + codes_dictionary[key] + '.fas ' + path_results + '\n')) 
        ## Print command to run process_bowtie_output.py to job file
        TASKSFORBC.write(str('perl ' + scriptpath + '/process_bowtie_output.pl ' + path_results + '\n'))
        if len(ptt)>0:
            TASKSFORBC.write(str('perl ' + scriptpath + '/normalize_processed.pl ' + path_results + '_processed.txt \n'))
            if arrayed == 1:
                if operon == 0:
                    TASKSFORBC.write(str('perl ' + scriptpath + '/map_genes_v2.py ' + path_indexes + '/' + str(ptt[0]) + ' ' + path_results + '_processed.txt_filter_cpm.txt ' + str(disruption) + '\n'))
                elif operon == 1:
                    TASKSFORBC.write(str('perl ' + scriptpath + '/map_genes_v1.pl ' + path_indexes + '/' + str(ptt[0]) + ' ' + path_results + '_processed.txt_filter_cpm.txt ' + str(disruption) + '\n'))
        TASKSFORBC.close()

        os.system(str('sbatch mappingjobs_' + myinput + '_' + key + '.job'))


elapsed = str(round((time.time() - start_time)/60))

print str("@ " + elapsed +  " minutes -- Jobs submitted!")        
            
        
                                  
    
