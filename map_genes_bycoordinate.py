#!/usr/bin/python -w -s


import sys
import glob


if len(sys.argv) < 3:
    print "\n path/map_genes_v2_py <.ptt> <.cpm file> <gene_disrupt_percent> (max=1) \n"

ptt = sys.argv[1]
file = sys.argv[2]
percent = float(sys.argv[3])
outfile = str(file + "_mapped")

#logfile = glob.glob('*scarf.log')
#LOG = open(logfile, 'w')

PTT = open(ptt, 'r')

dictionary = {}
for line in PTT:
    temp = line.strip().split('\t')
    if len(temp) > 1:
        sub = temp[0].strip().split('..')
        if len(sub) > 1:
            start = int(sub[0])
            end = int(sub[1])
            strand = temp[1]
            length = temp[2]
            gene = temp[5]
            synonym = temp[4]
            annotation = temp[8]
            if type(length) == str:
                length = end - start + 1
                dictionary[gene] = [start, end, strand, length, synonym, annotation]
            else:
                length = int(length * 3)
                dictionary[gene] = [start, end, strand, length, synonym, annotation]



PTT.close()

target = {}
for gene in dictionary.keys():
    start = int(dictionary[gene][0])
    end = int(dictionary[gene][1])
    length = int(dictionary[gene][3])
    newstart = round(start + (length*percent / int(2)),0)
    newend = round(end - (length*percent / int(2)),0)
    i = newstart - 1
    target[gene] = []
    while i >= newstart - 1 and i <= newend:
        target[gene].append(str(i).strip('.0'))
        i += 1

FILE = open(file, 'r')

hits = {}

for line in FILE:
    parsed = line.strip().split('\t')
    coordinate = parsed[1]
    count = float(parsed[2])
    for gene in target.keys():
        if coordinate in target[gene]:
            synonym = dictionary[gene][4]
            annotation = dictionary[gene][5]
            hits[coordinate] = [gene, synonym, annotation, count]               
                            

OUT = open(outfile, 'w')
#OUT.write(str('Gene ID\tCoordinate\tGene Name\tCount\tAnnotation\n'))
for coordinate in hits.keys():
    gene = str(hits[coordinate][0])
    synonym = str(hits[coordinate][1])
    annotation = str(hits[coordinate][2])
    count = str(hits[coordinate][3])
    coordinate = str(coordinate)
    OUT.write(str(gene + '\t' + coordinate + '\t' + synonym +'\t' + count + '\t' + annotation + '\n'))

OUT.close()
    

            
    
