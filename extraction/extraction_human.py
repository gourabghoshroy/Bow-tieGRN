#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np 

ids = {}
with open('..//datasets//human.node', 'r') as f:
    for line in f:
        linedata = line.rstrip().split("\t")
        val = linedata[1]
        if len(linedata) > 2:
            for i in range(2,len(linedata)):
                if len(linedata[i]) > len(val):
                    val = linedata[i]
        ids[linedata[0]] = val

edges = []
genes = []
ofile = open('..//networks//Human_gt.tsv', 'w')
with open('..//datasets//human.source', 'r') as f:
    for line in f:
        linedata = line.rstrip().split("\t")
        val1 = ids[linedata[1]]
        val2 = ids[linedata[3]]
        val = val1 + "\t" + val2
        if "MIR" in val1 or "hsa-let-" in val1 or "hsa-mir-" in val1 or "hsa-miR-" in val1:
            continue
        if val not in edges:
            edges.append(val)
            genes.append(val1)
            genes.append(val2)
            print(val,file=ofile)                

ofile.close()
genes = np.unique(genes)
print(len(genes))