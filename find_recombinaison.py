#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2020 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Search potential recombinaison from wgMLST database export"""

import sys
import os
import argparse
from lib import __version__
# import scipy as sc
# import scipy.stats as stats
# import matplotlib.backends.backend_pdf as pdf
# import matplotlib.pyplot as plt
# import math

desc = "Search potential gene recombinaison from wgMLST database export"
command = argparse.ArgumentParser(prog='find_recombinaison.py', \
    description=desc, usage='%(prog)s [options] genes alignment')
command.add_argument('-o', '--output', nargs='?', \
    type=argparse.FileType("w"), default=sys.stdout, \
    help='Output number of variations by genes (default:stdout)')
command.add_argument('genes', \
    type=argparse.FileType("r"), \
    help='genes list from mlst_extract_table with -e gene')
command.add_argument('alignment', \
    type=argparse.FileType("r"), \
    help='Multifasta alignment from mlst_export_sequence with -a')
command.add_argument('-v', '--version', action='version', version="pyMLST: "+__version__)


def compar_seqs(seqs):
    count = 0
    dim = len(seqs[0])
    for j in range(0, len(seqs[0])):
        d = set([s[j] for s in seqs])
        if '-' in d :
            d.remove('-')
        if len(d) > 1:
            count += 1
    return count

    
if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()    
    genes = [line.rstrip("\n") for line in args.genes]
    sys.stderr.write("Number of genes to look at : " + str(len(genes)) + "\n")
    output = args.output

    sequences = [[] for g in genes]
    samples = []

    ##load sequences by gene
    indice = 0
    for line in args.alignment:
        line = line.rstrip("\n")
        
        ##header
        if line.startswith(">"):
            indice = 0
            samples.append(line.lstrip(">"))
            continue

        ##check genes number correct
        if indice >= len(genes):
            raise Exception("The genes list seems not correspond to the alignment\n" + str(indice))
        
        ##genes
        sequences[indice].append(line)
        indice += 1

    ##check sequences are correctly align    
    for i,seqs in enumerate(sequences):
        if len(set([len(s) for s in seqs])) > 1:
            print(set([len(s) for s in seqs]))
            raise Exception("Following genes seems to be not align: " + genes[i])

    for i,seqs in enumerate(sequences):
        c = compar_seqs(seqs)
        output.write(genes[i] + "\t" + str(c) + "\t" + str(len(seqs[0])) + "\n")
        
