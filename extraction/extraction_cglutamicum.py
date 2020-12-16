#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np 

count = 0
genes = []
edges = []
with open('..//datasets//C_g_DSM_20300_=_ATCC_13032_regulations.csv', 'r') as f:
    for line in f:
        count += 1
        if count == 1:
            continue
        linedata = line.rstrip().split(",")
        val = linedata[2] + "\t" + linedata[6]
        if val not in edges:
            genes.append(linedata[2])
            genes.append(linedata[6])                
            edges.append(val)
                
genes = np.unique(genes)
print(len(genes))