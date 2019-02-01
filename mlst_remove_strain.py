#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Remove strains to an wgMLST database"""

import sys
import os
import argparse
import sqlite3
from lib import __version__

desc = "Remove strain to a wgMLST database and sequences specificaly associated"
command = argparse.ArgumentParser(prog='mlst_remove_strain.py', \
    description=desc, usage='%(prog)s [options] database strains')
command.add_argument('database', \
    type=argparse.FileType("r"), \
    help='Sqlite database of the wgMLST')
command.add_argument('strains', nargs='*',\
    type=str, help='List of strains to removed from database')
command.add_argument('-v', '--version', action='version', version="pyMLST: "+__version__)

if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()    
    database = args.database
    strains = args.strains

    if 'ref' in  args.strains:
        raise Exception("Ref schema could not be remove from this database")
    
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()
        
        for strain in strains:
            sys.stderr.write(strain + "     ")

            ## Search seq ids
            cursor.execute('''SELECT seqid FROM mlst WHERE souche=?''', (strain,))
            seqids = cursor.fetchall()
            if len(seqids) == 0:
                raise Exception("Strain name not found in database\n" + strain)
            
            ##remove sample
            cursor.execute('''DELETE FROM mlst WHERE souche=? ''', (strain,))
            
            ##remove seqs if no other strain have this seq
            for seqid in seqids:
                cursor.execute('''SELECT seqid FROM mlst WHERE seqid=?''', (seqid[0],))
                if cursor.fetchone() is None:
                    cursor.execute('''DELETE FROM sequences WHERE id=? ''', (seqid[0],))
            sys.stderr.write("OK\n")

        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
