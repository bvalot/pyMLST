#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL
import logging

import subprocess
import tempfile
from io import BytesIO

from pymlst.common.psl import Psl


class GenomePathNotFoundError(Exception):
    pass


def run_blat(path, genome, tmpfile, tmpout, identity, coverage):
    """Run Blat and return Psl Object"""
    command = [path, '-maxIntron=20', '-fine', '-minIdentity='+str(identity*100),
               genome.name, tmpfile.name, tmpout.name]
    proc = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    output, error = proc.communicate()
    for line in BytesIO(output).readlines():
        logging.info(line.decode().rstrip())
    have_error = False
    for line in BytesIO(error).readlines():
        have_error = True
        logging.error(line.decode().rstrip())
    if have_error:
        raise Exception('An error occurred while running BLAT,'
                        ' see the logs for more details.\n')
    genes = {}
    for line in open(tmpout.name, 'r'):
        try:
            int(line.split()[0])
        except (ValueError, IndexError):
            continue
        psl = Psl(line)
        if coverage <= psl.coverage <= 1:
            genes.setdefault(psl.gene_id(), []).append(psl)
    if len(genes) == 0:
        raise Exception("No path found for the coregenome")
    return genes


def blat_tmp():
    """Return a fasta and a psl temporary file"""
    tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)
    tmpout = tempfile.NamedTemporaryFile(mode='w+t', suffix='.psl', delete=False)
    tmpout.close()
    return tmpfile, tmpout
