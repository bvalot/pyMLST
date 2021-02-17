#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Create a wgMLST database"""

import sys

import time
import sqlite3
import os
import click
from Bio import SeqIO

from pymlst.lib import sql
from pymlst.wg_commands.db.database import Database


# def update_duplicate(cursor, gene):
#     cursor.execute('''SELECT id FROM sequences WHERE sequence=?''', (str(gene.seq).upper(),))
#     seqid = cursor.fetchone()[0]
#     sys.stderr.write("Concatenate gene " + gene.id + "\n")
#     cursor.execute('''SELECT gene FROM mlst WHERE seqid=?''', (seqid,))
#     othergene = cursor.fetchone()[0]
#     cursor.execute('''UPDATE mlst SET gene = ? WHERE seqid = ?''', (othergene+";"+gene.id, seqid))
#
#
# def remove_duplicate(cursor, toremove):
#     sys.stderr.write("Remove duplicate sequence :" + str(len(toremove)) + "\n")
#     for seq in toremove:
#         cursor.execute('''SELECT id FROM sequences WHERE sequence=?''', (seq,))
#         seqid = cursor.fetchone()[0]
#         cursor.execute('''DELETE FROM sequences WHERE id=? ''', (seqid,))
#         cursor.execute('''DELETE FROM mlst WHERE seqid=? ''', (seqid,))


@click.command(name="create_db")
@click.option('--concatenate', '-c',
              is_flag=True,
              help='Automatically concatenate genes with duplicated sequences')
@click.option('--remove', '-r',
              is_flag=True,
              help='Automatically remove genes with duplicated sequences')
@click.argument('coregene', type=click.File('r'))
@click.argument('database', type=click.File('w'))
def cli(coregene, database, concatenate, remove):
    """Create a wgMLST database from a template"""
    database.close()
    verybegin = time.time()
    db = Database(os.path.abspath(database.name), create=True)
    genes = set()
    toremove = set()
    try:
        #db = sqlite3.connect(database.name)
        # cursor = db.cursor()
        # cursor2 = db.cursor()
        # cursor.execute('''CREATE TABLE IF NOT EXISTS
        #                   sequences(id INTEGER PRIMARY KEY, sequence TEXT unique)''')
        # cursor.execute('''CREATE TABLE IF NOT EXISTS
        #                   mlst(id INTEGER PRIMARY KEY, souche TEXT, gene TEXT, seqid INTEGER)''')
        for gene in SeqIO.parse(coregene, 'fasta'):
            if gene.id in genes:
                raise Exception("Two sequences have the same gene ID : " + gene.id)
            else:
                genes.add(gene.id)
           # try:
                #seqid = sql.add_sequence(cursor, str(gene.seq))
            #first = time.time()
            seqid = db.add_sequence(str(gene.seq))
            #seq_t = time.time()

            if seqid == -1:
                if concatenate:
                    db.concatenate_gene(gene)
                elif remove:
                    toremove.add(str(gene.seq).upper())
                else:
                    raise Exception("Two genes have the same sequence " + gene.id +
                                    "\nUse -c or -r options to manage it")
            else:
                db.add_mlst(sql.ref, gene.id, seqid)
                #mlst_t = time.time()

                # print('Time seq: ' + str(seq_t - first))
                # print('Time mlst: ' + str(mlst_t - seq_t))

                # sql.add_mlst(cursor2, sql.ref, gene.id, seqid)
                # sql.add_mlst(sql.ref, gene.id, seqid)
            # except sqlite3.IntegrityError:  # Because sequence is 'unique'
            #     if concatenate:
            #         update_duplicate(cursor, gene)
            #     elif remove:
            #         toremove.add(str(gene.seq).upper())
            #     else:
            #         raise Exception("Two genes have the same sequence " + gene.id +
            #                         "\nUse -c or -r options to manage it")
        if toremove:
            db.remove_sequences(toremove)

        # Add index
        # sql.index_database(cursor)

        db.commit()

        veryend = time.time()

        print("Global time : ", str(veryend - verybegin))
    except Exception:
        #db.rollback()
        raise
    finally:
        db.close()
