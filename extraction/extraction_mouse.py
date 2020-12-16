#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np 

edges = []
genes = []

ofile = open('..//networks//Mouse_gt.tsv', 'w')
with open('..//datasets//mouse.source', 'r') as f:
    for line in f:
        linedata = line.rstrip().split("\t")
        val = linedata[0] + "\t" + linedata[2]
        if "mmu-let-" in linedata[0] or "mmu-miR-" in linedata[0]:
            continue
        if val not in edges:
            edges.append(val)
            genes.append(linedata[0])
            genes.append(linedata[2])
            print(val,file=ofile)

ofile.close()
genes = np.unique(genes)
print(len(genes))






