#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Remove gene to an wgMLST database"""

import sys
import os
import argparse
import sqlite3

desc = "Remove gene to an wgMLST database and sequences specificaly associated"
command = argparse.ArgumentParser(prog='mlst_remove_gene.py', \
    description=desc, usage='%(prog)s [options] database')
command.add_argument('-l', '--lists', \
    type=argparse.FileType("r"), \
    help='List of genes to removed on the wgMLST database')
command.add_argument('-g', '--gene', \
    type=str, \
    help='Gene(s) to removed on the wgMLST database. Multiples references could be add with space')
command.add_argument('database', \
    type=argparse.FileType("r"), \
    help='Sqlite database of the wgMLST')
    
if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()    
    database = args.database

    ##list genes to removed
    genes = []
    if args.lists is not None:
        for line in args.lists.readlines():
            genes.append(line.rstrip("\n"))
    if args.gene is not None:
        genes.extend(args.gene.split(" "))
    if len(genes) == 0:
        raise Exception("No gene to removed found.\n")
    genes = set(genes)
        
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()
        
        for gene in genes:
            sys.stderr.write(gene + "     ")

            ## Search seq ids
            cursor.execute('''SELECT seqid FROM mlst WHERE gene=?''', (gene,))
            seqids = cursor.fetchall()
            if len(seqids) == 0:
                raise Exception("Gene name not found in database\n" + gene)
            
            ##remove sample
            cursor.execute('''DELETE FROM mlst WHERE gene=? ''', (gene,))
            
            ##remove seqs if no other gene have this seq
            for seqid in seqids:
                cursor.execute('''SELECT seqid FROM mlst WHERE seqid=?''', (seqid[0],))
                if cursor.fetchone() is None:
                    cursor.execute('''DELETE FROM sequences WHERE id=? ''', (seqid[0],))
            sys.stderr.write("OK\n")

        db.commit()
    except Exception as e:
        db.rollback()
        raise e
    finally:
        db.close()
