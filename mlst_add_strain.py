#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Add a strain to wgMLST database"""

import sys
import os
import argparse
import sqlite3
from Bio import SeqIO
import shutil
import lib.psl as psl
import lib.blat as blat
from lib import __version__
import lib.sql as sql

desc = "Add a strain to the wgMLST database"
command = argparse.ArgumentParser(prog='mlst_add_strain.py', \
    description=desc, usage='%(prog)s [options] genome database')
command.add_argument('-s', '--strain', nargs='?', \
    type=str, default=None, \
    help='Name of the strain (default:genome name)')
command.add_argument('-i', '--identity', nargs='?', \
    type=float, default=0.95, \
    help='Minimun identity to search gene (default=0.95)')
command.add_argument('-c', '--coverage', nargs='?', \
    type=float, default=0.9, \
    help='Minimun coverage to search gene (default=0.9)')
command.add_argument('-p', '--path', nargs='?', \
    type=str, default="/usr/bin", \
    help='Path to BLAT executable (default=/usr/bin)')
command.add_argument('genome', \
    type=argparse.FileType("r"), \
    help='Genome of the strain')
command.add_argument('database', \
    type=argparse.FileType("r"), \
    help='Sqlite database to store MLST')
command.add_argument('-v', '--version', action='version', version="pyMLST: "+__version__)

def create_coregene(cursor, tmpfile):
    cursor.execute('''SELECT gene, sequence FROM mlst JOIN sequences 
                      ON mlst.seqid = sequences.id
                      WHERE mlst.souche = ?''', (sql.ref,))
    coregenes = []
    for row in cursor.fetchall():
        tmpfile.write('>' + row[0] + "\n" + row[1] + "\n")
        coregenes.append(row[0])
    return coregenes
    
def insert_sequence(cursor, sequence):
    try:
        return sql.add_sequence(cursor, sequence)
    except sqlite3.IntegrityError:
        cursor.execute('''SELECT id FROM sequences WHERE sequence=?''', (sequence.upper(),))
        return cursor.fetchone()[0]

def read_genome(genome):
    seqs = {}
    for seq in SeqIO.parse(genome, 'fasta'):
        seqs[seq.id] = seq
    return seqs
    
if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()    
    database = args.database
    genome = args.genome
    if args.identity<0 or args.identity > 1:
        raise Exception("Identity must be between 0 to 1")
    path = blat.test_blat_exe(args.path)

    name = args.strain
    if name is None:
        name = genome.name.split('/')[-1]
    if ";" in name:
        raise Exception("Strain name must not contains special ';'\n")

    tmpfile, tmpout = blat.blat_tmp()
    
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()
        cursor2 = db.cursor()

        ##index database for old one
        sql.index_database(cursor)

        ##verify strain not already on the database
        cursor.execute('''SELECT DISTINCT souche FROM mlst WHERE souche=?''', (name,))
        if cursor.fetchone() is not None:
            raise Exception("Strain is already present in database:\n"+name)
        
        ##read coregene
        coregenes = create_coregene(cursor, tmpfile)
        tmpfile.close()

        ##BLAT analysis
        sys.stderr.write("Search coregene with BLAT\n")
        genes = blat.run_blat(path, genome, tmpfile, tmpout, args.identity, args.coverage)
        sys.stderr.write("Finish run BLAT, found " + str(len(genes)) + " genes\n")
        
        ##add sequence MLST
        seqs = read_genome(genome)
        bad = 0
        for coregene in coregenes:
            if coregene not in genes:
                continue
            for gene in genes.get(coregene):
                seq = seqs.get(gene.chro, None)
                if seq is None:
                    raise Exception("Chromosome ID not found " + gene.chro)

                ##Correct coverage
                if gene.coverage != 1:
                    if gene.searchPartialCDS(seq, args.coverage) is False:
                        sys.stderr.write("Gene " + gene.geneId() + " partial: removed\n")
                        bad += 1
                        continue
                    else:
                        sys.stderr.write("Gene " + gene.geneId() + " fill: added\n")

                ##Verify CDS
                if psl.testCDS(gene.getSequence(seq), False) is False:
                    if gene.searchCorrectCDS(seq, args.coverage) is False:
                        sys.stderr.write("Gene " + gene.geneId() + " not correct: removed\n")
                        bad += 1
                        continue
                    else:
                        sys.stderr.write("Gene " + gene.geneId() + " correct: added\n")
 
                ##add sequence and MLST
                sequence = gene.getSequence(seq)
                
                ##Insert data in database
                seqid = insert_sequence(cursor, str(sequence))
                sql.add_mlst(cursor2, name, gene.geneId(), seqid)

        db.commit()
        sys.stderr.write("Add " + str(len(genes)-bad) + " new MLST gene to database\n")
        sys.stderr.write("FINISH\n")
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
        if os.path.exists(tmpfile.name):        
            os.remove(tmpfile.name)
        if os.path.exists(tmpout.name):        
            os.remove(tmpout.name)
