#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np 

count = 0
genes = []
edges = []
notvalid = ["NA","\N","N/A"]

ofile = open('..//networks//Arabidopsis_gt.tsv', 'w')
with open('..//datasets//AtRegNet.csv', 'r') as f:
    for line in f:
        count += 1
        if count == 1:
            continue
        linedata = line.rstrip().split(",")
        if len(linedata) < 8:
            continue
        for i in range(len(linedata)):
            linedata[i] = linedata[i].strip().upper()
        locindex = 4
        dirindex = 6
        if linedata[3] == "ENT1":
            locindex = 5
            dirindex = 7
        if linedata[dirindex][:6] != "DIRECT":
            continue
        if linedata[3] in notvalid or linedata[locindex] in notvalid:
            continue
        slocus = linedata[1].split("/")
        sname = linedata[0].split("/")
        for i in range(len(slocus)):
            if sname[i] in notvalid or slocus[i] in notvalid:
                continue
            val = slocus[i] + "\t" + linedata[locindex]
            if val not in edges:
                edges.append(val)
                genes.append(slocus[i])
                genes.append(linedata[locindex])
                print(val,file=ofile)               
                
ofile.close()               
genes = np.unique(genes)
print(len(genes))
            
            
