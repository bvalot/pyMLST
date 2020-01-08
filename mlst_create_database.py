#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Create a wgMLST database"""

import sys
import argparse
import sqlite3
from Bio import SeqIO
from lib import __version__
import lib.sql as sql

desc = "Create a wgMLST database from a template"
command = argparse.ArgumentParser(prog='mlst_create_database.py', \
    description=desc, usage='%(prog)s [options] coregene database')
command.add_argument('-c', '--concatenate', \
    action='store_true', \
    help='Automaticaly concatenate redondancy of genes with same sequence')
command.add_argument('-r', '--remove', \
    action='store_true', \
    help='Automaticaly remove genes with duplicate sequence')
command.add_argument('coregene', \
    type=argparse.FileType("r"), \
    help='Coregene fasta file as template of MLST')
command.add_argument('database', \
    type=argparse.FileType("w"), \
    help='Sqlite database to store MLST')
command.add_argument('-v', '--version', action='version', version="pyMLST: "+__version__)

def update_duplicate(cursor, gene):
    cursor.execute('''SELECT id FROM sequences WHERE sequence=?''', (str(gene.seq).upper(),))
    seqid = cursor.fetchone()[0]
    sys.stderr.write("Concatenate gene " + gene.id + "\n")
    cursor.execute('''SELECT gene FROM mlst WHERE seqid=?''', (seqid,))
    othergene = cursor.fetchone()[0]
    cursor.execute('''UPDATE mlst SET gene = ? WHERE seqid = ?''', (othergene+";"+gene.id, seqid))

def remove_duplicate(cursor, toremove):
    sys.stderr.write("Remove duplicate sequence :" + str(len(toremove)) + "\n")
    for seq in toremove:
        cursor.execute('''SELECT id FROM sequences WHERE sequence=?''', (seq,))
        seqid = cursor.fetchone()[0]
        cursor.execute('''DELETE FROM sequences WHERE id=? ''', (seqid,))
        cursor.execute('''DELETE FROM mlst WHERE seqid=? ''', (seqid,))  
    
if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()
    database = args.database
    database.close()
    genes = set()
    toremove = set()
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()
        cursor2 = db.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS
                          sequences(id INTEGER PRIMARY KEY, sequence TEXT unique)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS
                          mlst(id INTEGER PRIMARY KEY, souche TEXT, gene TEXT, seqid INTEGER)''')
        for gene in SeqIO.parse(args.coregene, 'fasta'):
            if gene.id in genes:
                raise Exception("Two sequences have the same gene ID : " + gene.id)
            else:
                genes.add(gene.id)
            try:
                seqid = sql.add_sequence(cursor, str(gene.seq))
                sql.add_mlst(cursor2, sql.ref, gene.id, seqid)
            except sqlite3.IntegrityError:
                if args.concatenate:
                    update_duplicate(cursor, gene)
                elif args.remove:
                    toremove.add(str(gene.seq).upper())  
                else:
                    raise Exception("Two gene have the same sequence " + gene.id + \
                                    "\nUse -c or -r options to manage it")
        if toremove:
            remove_duplicate(cursor, toremove)

        ## Add index
        sql.index_database(cursor)
        
        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
