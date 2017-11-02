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
print target
FILE = open(file, 'r')

hits = {}
for line in FILE:
    parsed = line.strip().split('\t')
    hit = parsed[1]
    count = float(parsed[2])
    unique = 1
    for gene in target.keys():
        if hit in target[gene]:
            synonym = dictionary[gene][4]
            annotation = dictionary[gene][5]
            if gene not in hits.keys():
                hits[gene] = [synonym, annotation, count, unique]               
            if gene in hits.keys():
                unique += 1
                hits[gene][2] += count
                hits[gene][3] = unique
print hits                
                

OUT = open(outfile, 'w')
OUT.write(str('Gene ID\tUnique insertions\tCount\tGene Name\t Annotation\n'))
for gene in hits.keys():
    gene = str(gene)
    synonym = str(hits[gene][0])
    annotation = str(hits[gene][1])
    count = str(hits[gene][2])
    unique = str(hits[gene][3])
    OUT.write(str(gene + '\t' + unique + '\t' + count + '\t' + synonym + '\t' + annotation + '\n'))

OUT.close()
    

            
    
