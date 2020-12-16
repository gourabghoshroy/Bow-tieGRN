#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np 

edges = []
genes = []
with open('..//datasets//trrust_rawdata.human.tsv', 'r') as f:
    for line in f:
        linedata = line.rstrip().split("\t")
        val = linedata[0] + "\t" + linedata[1]
        if val not in edges:
            edges.append(val)
            genes.append(linedata[0])
            genes.append(linedata[1])
                

genes = np.unique(genes)
print(len(genes))    