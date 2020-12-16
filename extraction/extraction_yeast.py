#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np 

genes = []
edges = []

ofile = open('..//networks//Yeast_gt.tsv', 'w')
with open('..//datasets//TRP_binding_network.txt', 'r') as f:
    for line in f:
        linedata = line.rstrip().split("\t")
        testval = linedata[3].split("->")
        if len(testval) == 1:
            val = testval[0][:len(testval[0])/2].upper() + "\t" + testval[0][len(testval[0])/2:].upper()
            if val not in edges:
                edges.append(val)
                genes.append(testval[0][:len(testval[0])/2].upper())
                genes.append(testval[0][len(testval[0])/2:].upper())
                print(val,file=ofile)
        for i in range(len(testval)-1):
            val = testval[i].upper() + "\t" + testval[i+1].upper()
            if val not in edges:
                edges.append(val)
                genes.append(testval[i].upper())
                genes.append(testval[i+1].upper())
                print(val,file=ofile)
                
ofile.close()               
genes = np.unique(genes)
print(len(genes))
            
            
