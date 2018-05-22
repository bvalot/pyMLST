#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Create a wgMLST database"""

import sys
import argparse
import sqlite3
from Bio import SeqIO

desc = "Create a wgMLST database from a template"
command = argparse.ArgumentParser(prog='mlst_create_strain.py', \
    description=desc, usage='%(prog)s [options] coregene database')
# command.add_argument('-s', '--strain', nargs='?', \
#     type=str, default=None, \
#     help='Name of the strain of coregene sequence (default:coregene name)')
command.add_argument('coregene', \
    type=argparse.FileType("r"), \
    help='Coregene fasta file as template of MLST')
command.add_argument('database', \
    type=argparse.FileType("w"), \
    help='Sqlite database to stock MLST')

if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()
    database = args.database
    # name = args.strain
    # if name is None:
    #     name = args.coregene.name.split('/')[-1]
    # print(name)
    name = "ref"
    database.close()
    genes = set()
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
                raise Exception("Two sequences have the same gene ID : " + seq.geneid)
            else:
                genes.add(gene.id)
            cursor.execute('''INSERT INTO sequences(sequence)
                              VALUES(?)''', (str(gene.seq),))
            cursor2.execute('''INSERT INTO mlst(souche, gene, seqid)
                              VALUES(?,?,?)''', (name, gene.id, cursor.lastrowid))
        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
