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

genes = []
with open('..//datasets//ENCODE', 'r') as f:
    for line in f:
        linedata = line.rstrip().split("\t")
        if linedata[0] == "79923":
            gene1 = "NANOG"
        elif linedata[0] not in ids.keys():
            continue
        else:
            gene1 = ids[linedata[0]]
        if "MIR" in gene1 or "hsa-let-" in gene1 or "hsa-mir-" in gene1 or "hsa-miR-" in gene1:
            continue
        genes.append(linedata[0])
        genes.append(linedata[1])
            
        
        
genes = np.unique(genes)
print(len(genes))