#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL
import logging
import sys
import subprocess
import tempfile
from pymlst.common.psl import Psl


def run_blat(path, genome, tmpfile, tmpout, identity, coverage):
    """Run Blat and return Psl Object"""
    command = [path, '-maxIntron=20', '-fine', '-minIdentity='+str(identity*100),\
               genome.name, tmpfile.name, tmpout.name]
    proc = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=sys.stderr)
    error = ""
    for line in proc.stderr:
        error += line.decode()
    if error != "":
        logging.error("An error occurred while running BLAT\n")
        raise Exception(error)
    genes = {}
    for line in open(tmpout.name, 'r'):
        try:
            int(line.split()[0])
        except:
            continue
        psl = Psl(line)
        if psl.coverage >=coverage and psl.coverage <= 1:
            genes.setdefault(psl.geneId(),[]).append(psl)
    if len(genes) == 0:
        raise Exception("No path found for the coregenome")
    return genes


def blat_tmp():
    """Return a fasta and a psl temporary file"""
    tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)
    tmpout = tempfile.NamedTemporaryFile(mode='w+t', suffix='.psl', delete=False)
    tmpout.close()
    return tmpfile, tmpout
