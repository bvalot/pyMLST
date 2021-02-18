#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Create a wgMLST database"""

import time
import os
import click
import sys
from Bio import SeqIO

from pymlst.lib import sql
from pymlst.wg_commands.db.database import Database


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
    to_remove = set()

    try:
        for gene in SeqIO.parse(coregene, 'fasta'):
            if gene.id in genes:
                raise Exception("Two sequences have the same gene ID: " + gene.id)
            else:
                genes.add(gene.id)

            added, seq_id = db.add_sequence(str(gene.seq))

            if not added:
                if concatenate:
                    db.concatenate_gene(seq_id, gene.id)
                    sys.stderr.write("Concatenate gene " + gene.id + "\n")
                elif remove:
                    to_remove.add(seq_id)
                else:
                    raise Exception("Two genes have the same sequence " + gene.id +
                                    "\nUse -c or -r options to manage it")
            else:
                db.add_mlst(sql.ref, gene.id, seq_id)

        if to_remove:
            db.remove_sequences(to_remove)
            sys.stderr.write("Remove duplicate sequence: " + str(len(to_remove)) + "\n")

        db.commit()

        veryend = time.time()

        print("Global time: ", str(veryend - verybegin))

    except Exception:
        db.rollback()
        raise

    finally:
        db.close()
