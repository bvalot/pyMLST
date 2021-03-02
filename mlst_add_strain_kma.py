#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Add a strain to wgMLST database from kma analysis"""

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

desc = "Add a strain to the wgMLST database from kma analysis"
command = argparse.ArgumentParser(prog='mlst_add_strain_kma.py', \
    description=desc, usage='%(prog)s [options] kma database')
command.add_argument('-s', '--strain', nargs='?', \
    type=str, default=None, \
    help='Name of the strain (default:genome name)')
command.add_argument('-i', '--identity', nargs='?', \
    type=float, default=0.95, \
    help='Minimun identity to search gene (default=0.95)')
command.add_argument('-c', '--coverage', nargs='?', \
    type=float, default=0.98, \
    help='Minimun alignement coverage to add gene (default=0.98)')
command.add_argument('-x', '--depth', nargs='?', \
    type=float, default=10, \
    help='Minimun reads depth to add gene (default=10)')
command.add_argument('kma', \
    type=argparse.FileType("r"), \
    help='KMA result file (.res)')
command.add_argument('database', \
    type=argparse.FileType("r"), \
    help='Sqlite database to store MLST')
command.add_argument('-v', '--version', action='version', version="pyMLST: "+__version__)

def get_all_gene(cursor):
    cursor.execute('''SELECT distinct(gene) FROM mlst WHERE souche = ?''', (sql.ref,))
    return [row[0] for row in cursor.fetchall()]

def insert_sequence(cursor, sequence):
    try:
        return sql.add_sequence(cursor, sequence)
    except sqlite3.IntegrityError:
        cursor.execute('''SELECT id FROM sequences WHERE sequence=?''', (sequence.upper(),))
        return cursor.fetchone()[0]

def read_fasta(fasta):
    seqs = {}
    for seq in SeqIO.parse(fasta, 'fasta'):
        seqs[seq.id] = seq.seq
    return seqs

def read_result(kma, cover, ident, reads):
    kmas = []
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
    
if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()    
    database = args.database
    kma = args.kma
    if args.identity<0 or args.identity > 1:
        raise Exception("Identity must be between 0 to 1")
    if args.coverage<0 or args.coverage > 1:
        raise Exception("Coverage must be between 0 to 1")
   # path = blat.test_blat_exe(args.path)

    name = args.strain
    if name is None:
        name = kma.name.split('/')[-1]
    if ";" in name:
        raise Exception("Strain name must not contains special ';'\n")

    fasta = kma.name.replace(".res", ".fsa")
    if os.path.exists(fasta) is False:
        raise Exception("Fsa file corresponding to kma results not found\n" + kma.name)
    #tmpfile, tmpout = blat.blat_tmp()
    
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()
        cursor2 = db.cursor()

        ##index database for old one
        sql.index_database(cursor)

        # ##verify strain not already on the database
        cursor.execute('''SELECT DISTINCT souche FROM mlst WHERE souche=?''', (name,))
        if cursor.fetchone() is not None:
            raise Exception("Strain is already present in database:\n"+name)
        
        # ##read coregene
        coregenes = get_all_gene(cursor)
        # tmpfile.close()

        ##read results
        kma_res = read_result(kma, args.coverage, args.identity, args.depth)
        seqs = read_fasta(fasta)
        
        ##add sequence MLST
        valid = 0
        minus = 0
        frame = 0
        duplicate = set() ##to prevent add sames alleles two times
        
        for res in kma_res:
            seq = seqs.get(res)
            if seq is None:
                raise Exception(res + " not found in the fasta files")

            ## test minus
            b = (seq.count('a') + seq.count('t') + seq.count('c') + \
                 seq.count('g'))
            if b != 0:
                minus +=1
                sys.stderr.write(res + " Remove incertain\n")
                continue

            ## test CDS
            try:
                seq.translate(cds=True, table=11)
            except:
                frame += 1
                sys.stderr.write(res + " Remove bad CDS\n")
                continue
 
            ##add sequence and MLST
            gene = res.split("_")[0]
            if gene not in coregenes:
                sys.stderr.write("WARNINGS: gene " + gene + " not present in database\n")
                continue
            valid +=1
            seqid = insert_sequence(cursor, str(seq))
            if str(seqid)+res.split("_")[0] not in duplicate:
                sql.add_mlst(cursor2, name, res.split("_")[0], seqid)
                duplicate.add(str(seqid)+res.split("_")[0])
            
        db.commit()
        sys.stderr.write("Add " + str(valid) + " new MLST genes to database\n")
        sys.stderr.write("Remove " + str(minus) + " genes with uncertain bases\n")
        sys.stderr.write("Remove " + str(frame) + " genes with bad CDS\n")  
        sys.stderr.write("FINISH\n")
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
