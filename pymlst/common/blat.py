#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL
import logging
import subprocess
import tempfile
import os
from io import BytesIO

from pymlst import config
from pymlst.common.psl import Psl
from pymlst.common import exceptions


def create_2bit_database(fasta_file, db_path):
    """Create a 2bit database from a fasta file for BLAT usage"""
    path_faToTwoBit = config.get_binary_path('faToTwoBit')
    if path_faToTwoBit is None:
        raise exceptions.BinaryNotFound('faToTwoBit binary was not found (part of UCSC tools)')
    
    two_bit_file = db_path + '.2bit'
    
    # Convert fasta to 2bit
    command = [path_faToTwoBit, fasta_file, two_bit_file]
    proc = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    output, error = proc.communicate()
    
    if proc.returncode != 0:
        raise exceptions.PyMLSTError(f'Error creating 2bit database: {error.decode()}')
    
    return two_bit_file


def run_blat_with_2bit(genome_path, two_bit_file, tmpout, identity, coverage, maxintron=20):
    """Run Blat using a 2bit database file and return Psl Object
    
    Args:
        genome_path: Path to the genome file
        two_bit_file: Path to the 2bit database file
        tmpout: Temporary output file for BLAT results
        identity: Minimum identity threshold
        coverage: Minimum coverage threshold
        maxintron: Maximum intron size
        
    Returns:
        Dictionary of genes with Psl objects
    """
    path = config.get_binary_path('blat')
    if path is None:
        raise exceptions.BinaryNotFound('BLAT binary was not found')

    # Run BLAT command
    command = [path, '-maxIntron='+str(maxintron), '-fine', 
               '-minIdentity='+str(identity*100), 
               '-noHead',  # Remove header to make parsing easier
               genome_path, two_bit_file, tmpout.name]
    
    logging.debug(f"Running BLAT command: {' '.join(command)}")
    
    proc = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    output, error = proc.communicate()
    
    for line in BytesIO(output).readlines():
        logging.debug(line.decode().rstrip())
    
    have_error = False
    for line in BytesIO(error).readlines():
        have_error = True
        logging.error(line.decode().rstrip())
    
    if have_error:
        raise exceptions.PyMLSTError('An error occurred while running BLAT')
    
    genes = {}
    for line in open(tmpout.name, 'r'):
        try:
            int(line.split()[0])
        except (ValueError, IndexError):
            continue
        psl_obj = Psl(line)
        if coverage <= psl_obj.coverage <= 1:
            genes.setdefault(psl_obj.gene_id(), []).append(psl_obj)
    
    if len(genes) == 0:
        raise exceptions.CoreGenomePathNotFound('No path was found for the core genome')
    
    return genes


def is_2bit_database(db_path):
    """Check if 2bit database exists"""
    two_bit_file = db_path + '.2bit'
    return os.path.exists(two_bit_file)


def delete_2bit_database(db_path):
    """Delete 2bit database files"""
    two_bit_file = db_path + '.2bit'
    if os.path.exists(two_bit_file):
        os.remove(two_bit_file)


def blat_tmp():
    """Return a fasta and a psl temporary file"""
    tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)
    tmpout = tempfile.NamedTemporaryFile(mode='w+t', suffix='.psl', delete=False)
    tmpout.close()
    return tmpfile, tmpout


# Keep original function for backward compatibility
def run_blat(genome, tmpfile, tmpout, identity, coverage, maxintron=20):
    """Original run_blat function for backward compatibility"""
    path = config.get_binary_path('blat')
    if path is None:
        raise exceptions.BinaryNotFound('BLAT binary was not found')

    command = [path, '-maxIntron='+str(maxintron), '-fine', \
               '-minIdentity='+str(identity*100), \
               genome.name, tmpfile.name, tmpout.name]
    proc = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

    output, error = proc.communicate()
    for line in BytesIO(output).readlines():
        logging.debug(line.decode().rstrip())
    have_error = False
    for line in BytesIO(error).readlines():
        have_error = True
        logging.error(line.decode().rstrip())
    if have_error:
        raise exceptions.PyMLSTError(
            'An error occurred while running BLAT')
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
        raise exceptions.CoreGenomePathNotFound(
            'No path was found for the core genome')
    return genes
