#!/usr/bin/python -w -s

import sys
import glob

filepath = sys.argv[1]
allmapped = glob.glob(str(sys.argv[1] + '/*mapped'))
filenumber = len(allmapped)

#coordinates = {}
name = {}

for file in allmapped:
    IN = open(file, 'r')
    for line in IN:
	unique = 1
        temp = line.strip().split('\t')
#        bp = temp[2]
        synonym = temp[2]
        gene = temp[0]
        annotation = temp[4]
        
#        if bp not in coordinates.keys():
#            coordinates[bp] = [gene, synonym, annotation]
        if synonym not in name.keys():
            name[synonym] = [gene, synonym, annotation, unique]
	if synonym in name.keys():
	    name[synonym][3] += 1
    IN.close()

for file in allmapped:
    IN = open(file,'r')
    outfile = str(file + '_all.txt')
    OUT = open(outfile, 'w')
    OUT.write(str('Gene ID\Unique Insertions\tGene Name\tCount\tAnnotation\n'))
    allin = name.keys()
    print len(allin)
    for line in IN:
        temp = line.strip().split('\t')
#        coordinate = str(temp[1])
        synonym = str(temp[2])
 	gene = str(name[synonym][0])
	unique = str(name[synonym][3])
        count = str(temp[3])
        annotation = str(name[synonym][2])
        if synonym in allin:
            OUT.write(str(gene + '\t' + unique + '\t' + synonym +'\t' + count + '\t' + annotation + '\n'))
            allin.remove(synonym)
    for item in allin:
        gene = str(name[item][0])
        synonym = str(name[item][1])
        count = str(0)
        annotation = str(name[item][2])
        unique = str(name[item][3])
        OUT.write(str(gene + '\t' + unique + '\t' + synonym +'\t' + count + '\t' + annotation + '\n'))


    print len(allin)
    OUT.close()
        
        
        



