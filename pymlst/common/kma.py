#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2021 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL
import logging

import subprocess
import tempfile
from io import BytesIO
import os

from pymlst import config
from pymlst.common import utils
from pymlst.common.psl import Psl
from pymlst.common import exceptions

index = [".comp.b", ".length.b", ".name", ".seq.b"]
suffix = ".kma"

def run_kma(fastq, basename, identity, coverage, reads):
    """Run kma on fastq(s) and return sequences"""
    if is_database_indexing(basename) is False:
        raise exceptions.PyMLSTError('Dabatase must be index with KMA')
    
    path = config.get_binary_path('kma')
    if path is None:
        raise exceptions.BinaryNotFound('KMA binary was not found')

    with tempfile.NamedTemporaryFile('w+t') as tmp:
        baseout = tmp.name
    command = [path, '-t_db', basename+suffix, '-o', baseout, '-nf']
    if len(fastq) == 1:
        command.extend(['-i', fastq[0].name])
    elif len(fastq) == 2:
        command.extend(['-ipe', fastq[0].name, fastq[1].name])
    else:
        raise exceptions.PyMLSTError('Too many fastq files in input of run_kma')

    logging.info("Running KMA with cg/wgMLST database")
    proc = subprocess.Popen(command, stderr=subprocess.PIPE, \
                            stdout=subprocess.PIPE)

    output, error = proc.communicate()
    if os.path.exists(baseout + ".res") and os.path.exists(baseout + ".fsa"):
        for line in BytesIO(error).readlines():
            logging.debug(line.decode().rstrip())
    else:
        for line in BytesIO(error).readlines():
            logging.error(line.decode().rstrip())
        raise exceptions.PyMLSTError(
            'An error occurred while running KMA')   

    with open(baseout + ".res", 'r') as kma:
        kma_res = read_kma_res(kma, coverage, identity, reads)
    seqs = utils.read_genome(baseout + ".fsa")

    del_kma_tmp(baseout)
    if len(kma_res) == 0:
        raise exceptions.CoreGenomePathNotFound(
            'No path was found for the core genome')
    return kma_res,seqs


def del_kma_tmp(baseout):
    """Delete temporary file create by kma"""
    for a in [".aln", ".res", ".fsa"]:
        if os.path.exists(baseout + a):
            os.remove(baseout + a)
    
def is_database_indexing(basename):
    """Verify if a pyMLST database is indexing"""
    for i in index:
        if os.path.exists(basename + suffix + i) is False:
            return False
    return True

def index_database(basename, coregenes):
    """Index a database with kma if the base is not already indexing
    
    :coregene is a temporary file containing coregenes sequences
    """
    if is_database_indexing(basename) is False:
        path = config.get_binary_path('kma')
        if path is None:
            raise exceptions.BinaryNotFound('KMA binary was not found')
        logging.info("Indexing database %s with kma", \
                     os.path.basename(basename))
        
        command = [path, 'index', '-i', coregenes.name, '-o', basename + suffix]
        proc = subprocess.Popen(command, stderr=subprocess.PIPE, \
                                stdout=subprocess.PIPE)
        output, error = proc.communicate()
        if is_database_indexing(basename) is False:
            for line in BytesIO(error).readlines():
                logging.error(line.decode().rstrip())
            raise exceptions.PyMLSTError(
                'An error occurred while indexing KMA')
        else:
            for line in BytesIO(error).readlines():
                logging.debug(line.decode().rstrip())


def delete_indexing(basename):
    """Remove indexing file"""
    for i in index:
        if os.path.exists(basename + suffix + i):
            os.remove(basename + suffix + i)

def read_kma_res(kma, cover, ident, reads):
    kmas=[]
    header = kma.readline().rstrip("\n").split("\t")
    if len(header) != 11 or header[0].startswith("#Template") is False:
        raise Exception(kma.name + " seems not to be a kma result file\n")
    for line in kma:
        values = line.rstrip("\n").split("\t")
        if len(values) != 11:
            raise Exception("Incorrect line\n" + line)
        ele = {a:b.strip(" ") for a,b in zip(header, values)}
        if float(ele.get("Template_Coverage")) >= cover*100 and \
           float(ele.get("Query_Coverage")) >= cover*100 and \
           float(ele.get("Template_Identity")) >= ident*100 and \
           float(ele.get("Query_Identity")) >= ident*100 and \
           float(ele.get("Depth")) >= reads :
            kmas.append(ele.get("#Template"))
    return kmas

def index_tmpfile():
    return tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta')
