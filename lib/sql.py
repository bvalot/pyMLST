#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

import sqlite3

ref = "ref"

def index_database(cursor):
    cursor.execute('''CREATE INDEX IF NOT EXISTS ID_souche ON mlst (souche)''')
    cursor.execute('''CREATE INDEX IF NOT EXISTS ID_gene ON mlst (gene)''')
    cursor.execute('''CREATE INDEX IF NOT EXISTS ID_seqid ON mlst (seqid)''')

def add_sequence(cursor, seq):
    """Add sequence to database and return seqid""" 
    cursor.execute('''INSERT INTO sequences(sequence)
                      VALUES(?)''', (seq.upper(),))
    return cursor.lastrowid

def add_mlst(cursor, souche, gene, seqid):
    """Add MLST information to database"""
    cursor.execute('''INSERT INTO mlst(souche, gene, seqid)
                       VALUES(?,?,?)''', (souche, gene, seqid))
