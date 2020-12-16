#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np

genes = []
edges = []

ofile = open('..//networks//Ecoli_gt.tsv', 'w')
with open('..//datasets//network_tf_gene.txt', 'r') as f:
    for line in f:
        if line[0] == "#":
            continue
        lineData = line.rstrip().split("\t")
        if lineData[0] == "H-NS":
            lineData[0] = "hns"
        elif lineData[0] == "IHF":
            val = "ihfA\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("ihfA")
                genes.append(lineData[1])
                print(val,file=ofile) 
            val = "ihfB\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("ihfB")
                genes.append(lineData[1])
                print(val,file=ofile) 
            continue
        elif lineData[0] == "HU":
            val = "hupA\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("hupA")
                genes.append(lineData[1])
                print(val,file=ofile) 
            val = "hupB\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("hupB")
                genes.append(lineData[1])
                print(val,file=ofile) 
            continue
        elif lineData[0] == "FlhDC":
            val = "flhD\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("flhD")
                genes.append(lineData[1])
                print(val,file=ofile) 
            val = "flhC\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("flhC")
                genes.append(lineData[1])
                print(val,file=ofile) 
            continue
        elif lineData[0] == "HipAB":
            val = "hipA\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("hipA")
                genes.append(lineData[1])
                print(val,file=ofile) 
            val = "hipB\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("hipB")
                genes.append(lineData[1])
                print(val,file=ofile) 
            continue
        elif lineData[0] == "HigBA":
            val = "higA\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("higA")
                genes.append(lineData[1])
                print(val,file=ofile) 
            val = "higB\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("higB")
                genes.append(lineData[1])
                print(val,file=ofile) 
            continue
        elif lineData[0] == "RcsAB":
            val = "rcsA\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("rcsA")
                genes.append(lineData[1])
                print(val,file=ofile) 
            val = "rcsB\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append("rcsB")
                genes.append(lineData[1])
                print(val,file=ofile) 
            continue            
        elif "-" in lineData[0]:
            val2 = lineData[0].split("-")
            val2[0] = val2[0][0].lower()+val2[0][1].lower()+val2[0][2].lower()+val2[0][3:]
            val = val2[0] + "\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append(val2[0])
                genes.append(lineData[1])
                print(val,file=ofile)
            val2[1] = val2[1][0].lower()+val2[1][1].lower()+val2[1][2].lower()+val2[1][3:]
            val = val2[1] + "\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append(val2[1])
                genes.append(lineData[1])
                print(val,file=ofile)
            continue
        else:
            lineData[0] = lineData[0][0].lower()+lineData[0][1].lower()+lineData[0][2].lower()+lineData[0][3:]
        val = lineData[0] + "\t" + lineData[1]
        if val not in edges:
            edges.append(val)
            genes.append(lineData[0])
            genes.append(lineData[1])
            print(val,file=ofile)
                

with open('..//datasets//network_sigma_gene.txt', 'r') as f:
    for line in f:
        if line[0] == "#":
            continue
        lineData = line.rstrip().split("\t")
        sigfactors = lineData[0].split(",")
        for sigf in sigfactors:
            sigf = sigf.strip()
            if sigf == "Sigma24":
                sigf = "rpoE"
            elif sigf == "Sigma28":
                sigf = "fliA"
            elif sigf == "Sigma32":
                sigf = "rpoH"
            elif sigf == "Sigma70":
                sigf = "rpoD"
            elif sigf == "Sigma38":
                sigf = "rpoS"
            elif sigf == "Sigma54":
                sigf = "rpoN"
            val = sigf + "\t" + lineData[1]
            if val not in edges:
                edges.append(val)
                genes.append(sigf)
                genes.append(lineData[1])
                print(val,file=ofile)
                

ofile.close()       
genes = np.unique(genes)
print(len(genes))