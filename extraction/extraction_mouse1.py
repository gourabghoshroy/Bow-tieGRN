#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np 

count = 0
edges = []
genes = []
with open('..//datasets//ORTI-Sept-2016.tsv', 'r') as f:
    for line in f:
        count += 1
        if count == 1:
            continue
        linedata = line.rstrip().split("\t")
        if linedata[9].strip() != "Mus musculus":
            continue
        val = linedata[2] + "\t" + linedata[6]
        if linedata[2] == "-" or linedata[6] == "-":
            continue
        if val not in edges:
            edges.append(val)
            genes.append(linedata[2])
            genes.append(linedata[6])
                

genes = np.unique(genes)
print(len(genes))       