#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2020 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Search group of strains at a distance treshold"""

import sys
import os
import argparse
from lib import __version__
import networkx as nx
#import matplotlib.pyplot as plt
#import itertools as it

desc = "Search group os strain at a distance treshold"
command = argparse.ArgumentParser(prog='find_subgraph.py', \
    description=desc, usage='%(prog)s [options] distance')
command.add_argument('-o', '--output', nargs='?', \
    type=argparse.FileType("w"), default=sys.stdout, \
    help='Output result files (default:stdout)')
command.add_argument('-t', '--threshold', nargs='?', \
    type=int, default=50, \
    help='Minimum distance to conserve for extraction of group (default:50)')
command.add_argument('-e', '--export', nargs='?', \
    choices=['group', 'list', 'count'], default="group", \
    help='Defined the type of export (default:group)')
command.add_argument('distance', \
    type=argparse.FileType("r"), \
    help='Distance matrix from mlst_export_table with -e distance')
command.add_argument('-v', '--version',
    action='version', version="pyMLST: "+__version__)
    

if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()    
    output = args.output
    export = args.export
    threshold = args.threshold
    distance = args.distance

    # treshold=200
    # distance=open("distance-pyocls-m257", 'r')
    # group=open("distance-pyocls-m257.groups", 'w')
    # count=open("distance-pyocls-m257.count", 'w')

    ##load data
    samps = []
    dists = []
    try:
        strains = int(distance.readline().rstrip("\n"))
    except:
        raise Exception("The distance file seems not correctly formatted\n Not integer on first line")

    for line in distance.readlines():
        h = line.rstrip("\n").split("\t")
        samps.append(h[0])
        dists.append(h[1:])
        
    if len(samps) != strains:
        raise Exception("The distance file seems not correctly formatted\n Number of strains "+ str(len(samps)) + " doesn't correspond to " + str(strains))

    ##create graph
    G = nx.Graph()
    G.add_nodes_from(samps)

    for i,s in enumerate(samps):
        for j,d in enumerate(dists[i]):
            d = int(d)
            if i==j or d > threshold:
                continue
            G.add_edge(samps[i], samps[j], weight=d)

    ##extract interconnected subgraph
    ##count sample not found
    samps2 = set(samps)
    grps = []
    for subG in [G.subgraph(c) for c in nx.connected_components(G)]:

        inds = []
        for n in subG.nodes():
            samps2.remove(n)
            inds.append(samps.index(n))
        grps.append(inds)
    
    grps.sort(key=len,reverse=True)

    ##write result
    if export == "group":
        for i,g in enumerate(grps):
            #a = len(samps)*[0]
            output.write("Group" + str(i))
            for n in g:
                #a[n] = 1
                output.write(" " + samps[n])
            output.write("\n")
            
    elif export == "count":
        output.write("Group\t" + "\t".join(samps) + "\n")
        for i,g in enumerate(grps):
            a = len(samps)*[0]
            #group.write("Group" + str(i))
            for n in g:
                a[n] = 1
                #group.write(" " + samps[n])
            output.write(str(i) + "\t" + "\t".join(map(str,a))+ "\n")
            #group.write("\n")
    else:
       for i,g in enumerate(grps):
            #a = len(samps)*[0]
            #group.write("Group" + str(i))
            for n in g:
                #a[n] = 1
                output.write("Group" + str(i) + "\t" + samps[n] + "\n")
                #group.write(" " + samps[n])
            #write_count(count, str(i) + "\t" + "\t".join(map(str,a))+ "\n")
            #group.write("\n")
