#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import networkx as nx

#np.random.seed(100)
species = "Ecoli"    #Species names should match the name of the GRN file
randnum = 1000       #Number of networks with edge addition
numthr = 0.1         #Fraction of original edges to be added
#Output file to write the results - Number of nodes in each of the 7 bow-tie layers followed by number of regulators in each of the 7 layers 
outfilename = '..//output//'+species+'_bowtierandom'+str(int(numthr*100))+'addedges.tsv'  

genes = []
edges = []
regs = []      
with open('..//networks//'+species+'_gt.tsv', 'r') as f:   #Read GRN 
    for line in f:
        lineData = line.rstrip().split("\t")
        edges.append(lineData[0]+"\t"+lineData[1])
        genes.append(lineData[0])
        genes.append(lineData[1])
        regs.append(lineData[0])

genes = np.unique(genes)
regs = set(regs)

#Create directed graph 
net = nx.DiGraph()

for gene in genes:
    net.add_node(gene)
    
for edge in edges:
    lineData = edge.split("\t")
    net.add_edge(lineData[0], lineData[1])

ethr = np.int(np.round(numthr*len(edges)))

ofile = open(outfilename, 'w')
for i in range(randnum):     #For each network with random edge addition
    netn = net.copy()
    ecount = 0
    while (1>0):
        [gene1] = np.random.choice(list(regs),1)
        [gene2] = np.random.choice(genes,1)
        if not netn.has_edge(gene1,gene2):
            netn.add_edge(gene1,gene2)
            ecount += 1
            if ecount >= ethr:
                break
    
    # Bow-tie decomposition of the new network 
    count = 0
    lsc2 = set()      
    for c in sorted(nx.strongly_connected_components(netn),key=len,reverse=True):
        count += 1
        if count == 1:
            core = c
        elif count == 2:
            lsc2 = c
        else:
            break

    for node in core:
        dfsd = nx.dfs_successors(netn.reverse(),node)
        break

    dfsnodesgt = []
    for nodes in dfsd.values():
        dfsnodesgt =  dfsnodesgt + nodes
    dfsnodesgt = set(dfsnodesgt)
    inp = dfsnodesgt - core

    for node in core:
        dfsd = nx.dfs_successors(netn,node)
        break

    dfsnodesg = []
    for nodes in dfsd.values():
        dfsnodesg = dfsnodesg + nodes
    dfsnodesg = set(dfsnodesg)
    out = dfsnodesg - core

    intendrils  = set()
    outtendrils  = set()
    tubes = set()
    others = set()

    for node in netn.nodes():
        if node in core or node in inp or node in out:
            continue
        dfsd = nx.dfs_successors(netn.reverse(),node)
        dfsnodesgt = []
        for nodes in dfsd.values():
            dfsnodesgt = dfsnodesgt + nodes
        dfsnodesgt = set(dfsnodesgt)
        irv = False
        irvval = inp.intersection(dfsnodesgt)
        if len(irvval) != 0:
            irv = True    
        dfsd = nx.dfs_successors(netn,node)
        dfsnodes = []
        for nodes in dfsd.values():
            dfsnodes = dfsnodes + nodes
        dfsnodes = set(dfsnodes)
        vro = False
        vroval = out.intersection(dfsnodes)
        if len(vroval) != 0:
            vro = True
    
        if irv and not vro:
            intendrils.add(node)
        elif not irv and vro:
            outtendrils.add(node)
        elif irv and vro:
            tubes.add(node)
        else:
            others.add(node)
                
    corereg = core.intersection(regs)
    inpreg = inp.intersection(regs)
    outreg = out.intersection(regs)
    intendrilsreg = intendrils.intersection(regs)
    outtendrilsreg = outtendrils.intersection(regs)
    tubesreg = tubes.intersection(regs)
    othersreg = others.intersection(regs)
                
    print(str(len(core))+"\t"+str(len(inp))+"\t"+str(len(out))+"\t"+str(len(intendrils))+"\t"+str(len(outtendrils))+"\t"+str(len(tubes))+"\t"+str(len(others))+"\t",end='',file=ofile)
    print(str(len(corereg))+"\t"+str(len(inpreg))+"\t"+str(len(outreg))+"\t"+str(len(intendrilsreg))+"\t"+str(len(outtendrilsreg))+"\t"+str(len(tubesreg))+"\t"+str(len(othersreg)),file=ofile)
                
ofile.close()
                

