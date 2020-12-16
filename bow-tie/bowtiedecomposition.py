#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import networkx as nx

species = "Ecoli"   #Species names should match the name of the GRN file

genes = []
edges = []
regs = []      
with open('..//networks//'+species+'_gt.tsv', 'r') as f: #Read GRN 
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
    
count = 0
lsc2 = set()       #2nd LSC
for c in sorted(nx.strongly_connected_components(net),key=len,reverse=True):
    count += 1
    if count == 1:
        core = c     #CORE LAYER
    elif count == 2:
        lsc2 = c
    else:
        break

for node in core:
    dfsd = nx.dfs_successors(net.reverse(),node)
    break

dfsnodesgt = []
for nodes in dfsd.values():
    dfsnodesgt =  dfsnodesgt + nodes
dfsnodesgt = set(dfsnodesgt)
inp = dfsnodesgt - core      #IN LAYER

for node in core:
    dfsd = nx.dfs_successors(net,node)
    break

dfsnodesg = []
for nodes in dfsd.values():
    dfsnodesg = dfsnodesg + nodes
dfsnodesg = set(dfsnodesg)
out = dfsnodesg - core       #OUT LAYER

intendrils  = set()
outtendrils  = set()
tubes = set()
others = set()

for node in net.nodes():     #REMAINING NODES IN TENDRILS, TUBES OR OTHERS LAYERS
    if node in core or node in inp or node in out:
        continue
    dfsd = nx.dfs_successors(net.reverse(),node)
    dfsnodesgt = []
    for nodes in dfsd.values():
        dfsnodesgt = dfsnodesgt + nodes
    dfsnodesgt = set(dfsnodesgt)
    irv = False
    irvval = inp.intersection(dfsnodesgt)
    if len(irvval) != 0:
        irv = True    
    dfsd = nx.dfs_successors(net,node)
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


#PERCENTAGE OF ALL NODES
print("Size wrt all nodes")
print("CORE size : "+str(len(core)))
print("CORE size % : "+str(len(core)*100.0/len(genes)))
print("2nd LSC size : "+str(len(lsc2)))
print("IN size : "+str(len(inp)))
print("IN size % : "+str(len(inp)*100.0/len(genes)))
print("OUT size : "+str(len(out)))
print("OUT size % : "+str(len(out)*100.0/len(genes)))    
print("INTENDRILS size : "+str(len(intendrils)))
print("INTENDRILS size % : "+str(len(intendrils)*100.0/len(genes)))
print("OUTTENDRILS size : "+str(len(outtendrils)))
print("OUTTENDRILS size % : "+str(len(outtendrils)*100.0/len(genes)))
print("TUBES size : "+str(len(tubes)))
print("TUBES size % : "+str(len(tubes)*100.0/len(genes)))
print("OTHERS size : "+str(len(others)))
print("OTHERS size % : "+str(len(others)*100.0/len(genes)))

print("-------")
#PERCENTAGE OF ALL REGULATORS
corereg = core.intersection(regs)
inpreg = inp.intersection(regs)
outreg = out.intersection(regs)
intendrilsreg = intendrils.intersection(regs)
outtendrilsreg = outtendrils.intersection(regs)
tubesreg = tubes.intersection(regs)
othersreg = others.intersection(regs)
print("Size wrt all regulators")
print("CORE size : "+str(len(corereg)))
print("CORE size % : "+str(len(corereg)*100.0/len(regs)))
print("IN size : "+str(len(inpreg)))
print("IN size % : "+str(len(inpreg)*100.0/len(regs)))
print("OUT size : "+str(len(outreg)))
print("OUT size % : "+str(len(outreg)*100.0/len(regs)))    
print("INTENDRILS size : "+str(len(intendrilsreg)))
print("INTENDRILS size % : "+str(len(intendrilsreg)*100.0/len(regs)))
print("OUTTENDRILS size : "+str(len(outtendrilsreg)))
print("OUTTENDRILS size % : "+str(len(outtendrilsreg)*100.0/len(regs)))
print("TUBES size : "+str(len(tubesreg)))
print("TUBES size % : "+str(len(tubesreg)*100.0/len(regs)))
print("OTHERS size : "+str(len(othersreg)))
print("OTHERS size % : "+str(len(othersreg)*100.0/len(regs)))


            
