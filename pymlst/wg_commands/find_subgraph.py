#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2020 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Search group of strains at a distance treshold"""

import sys
import click
import networkx as nx


desc = "Search group os strain at a distance treshold"


def write_count(count, texte):
    if count:
        count.write(texte)


@click.command()
@click.option('--output', '-o',
              type=click.File('w'),
              default=sys.stdout,
              help='Output group files (default:stdout)')
@click.option('--threshold', '-t',
              type=int, default=50,
              help='Minimum distance to conserve '
                   'for extraction of group (default:50)')
@click.option('--count', '-c',
              type=click.File('w'),
              help='Output count file')
@click.argument('distance',
                type=click.File('r'))
def cli(output, threshold, count, distance):
    """Search group os strain at a distance treshold"""

    # treshold=200
    # distance=open("distance-pyocls-m257", 'r')
    # group=open("distance-pyocls-m257.groups", 'w')
    # count=open("distance-pyocls-m257.count", 'w')

    # load data
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

    # create graph
    G = nx.Graph()
    G.add_nodes_from(samps)

    for i,s in enumerate(samps):
        for j,d in enumerate(dists[i]):
            d = int(d)
            if i==j or d > threshold:
                continue
            G.add_edge(samps[i], samps[j], weight=d)

    # extract interconnected subgraph
    # count sample not found
    samps2 = set(samps)
    grps = []
    for subG in [G.subgraph(c) for c in nx.connected_components(G)]:

        inds = []
        for n in subG.nodes():
            samps2.remove(n)
            inds.append(samps.index(n))
        grps.append(inds)
    
    grps.sort(key=len,reverse=True)

    # write result
    write_count(count, "Group\t" + "\t".join(samps) + "\n")
    for i,g in enumerate(grps):
        a = len(samps)*[0]
        output.write("Group" + str(i))
        for n in g:
            a[n] = 1
            output.write(" " + samps[n])
        write_count(count, str(i) + "\t" + "\t".join(map(str,a))+ "\n")
        output.write("\n")

    if count:
        count.close()
    output.close()
