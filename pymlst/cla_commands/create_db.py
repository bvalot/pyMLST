#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Create a classical MLST database"""

import sys
import click
import sqlite3
from Bio import SeqIO

desc = "Create a classical MLST database from a sheme"


@click.command()
@click.argument('database',
                type=click.File('w'))
@click.argument('scheme',
                type=click.File('r'))
@click.argument('alleles',
                type=click.File('r'), nargs=-1, required=True)
def cli(database, scheme, alleles):
    """Create a classical MLST database from a sheme"""
    name = "ref"
    database.close()
    
    # Verify sheme list with fasta files
    header = scheme.readline().rstrip("\n").split("\t")
    if len(header) != len(alleles)+1:
        raise Exception("The number of genes in sheme don't correspond to the number of fasta file\n" \
                        + " ".join(header)  + "\n")
    fastas = {}
    for f in alleles:
        name = f.name.split("/")[-1]
        name = name[:name.rfind(".")]
        if name not in header:
            raise Exception("Gene " + name + " not found in sheme\n"+ " ".join(header))
        fastas[name]=f
        
    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()
        cursor2 = db.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS
                      sequences(id INTEGER PRIMARY KEY, sequence TEXT unique, gene TEXT, allele INTEGER)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS
                          mlst(id INTEGER PRIMARY KEY, st INTEGER, gene TEXT, allele INTEGER)''')
        # load sequence allele
        alleles = {}
        for g,f in fastas.items():
            alleles[g] = set()
            for seq in SeqIO.parse(f, 'fasta'):
                try:
                    if len(seq.id.split("_"))==2:
                        allele = int(seq.id.split("_")[1])
                    elif len(seq.id.split("-"))==2:
                        allele = int(seq.id.split("-")[1])
                    elif f.name.split(".")[0] in seq.id:
                        allele = int(seq.id.replace(f.name.split(".")[0], ""))
                    else:
                        allele = int(seq.id)
                except Exception:
                    raise Exception("Unable to obtain allele number for the sequence: " + seq.id)
                cursor.execute('''INSERT INTO sequences(sequence,gene,allele)
                              VALUES(?,?,?)''', (str(seq.seq).upper(),g,allele))
                alleles.get(g).add(allele)
        
        # load MLST sheme
        for line in scheme:
            h = line.rstrip("\n").split("\t")
            st = int(h[0])
            for g,a in zip(header[1:],h[1:]):
                if int(a) not in alleles.get(g):
                    sys.stderr.write("Unable to find the allele number "+a+" for gene "+g+ "; replace by 0\n")
                    # raise Exception("Unable to find the allele number "+a+" for gene "+g+ "\n" )
                    cursor2.execute('''INSERT INTO mlst(st, gene, allele)
                              VALUES(?,?,?)''', (st, g, 0))
                else:
                    cursor2.execute('''INSERT INTO mlst(st, gene, allele)
                              VALUES(?,?,?)''', (st, g, int(a)))
        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
