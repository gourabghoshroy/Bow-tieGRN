#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np 

count = 0
genes = []
edges = []

ofile = open('..//networks//Drosophila_gt.tsv', 'w')
with open('..//datasets//tf_gene.txt', 'r') as f:
    for line in f:
        count += 1
        if count == 1:
            continue
        linedata = line.rstrip().split("\t")
        val = linedata[0] + "\t" + linedata[1]
        if val not in edges:
            edges.append(val)
            genes.append(linedata[0])
            genes.append(linedata[1])
            print(val,file=ofile) 
                
ofile.close()                
genes = np.unique(genes)
print(len(genes))