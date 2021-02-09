#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Benoit Valot
# benoit.valot@univ-fcomte.fr
# UMR 6249 Chrono-Environnement, Besançon, France
# Licence GPL

"""Remove gene to a wgMLST database"""

import sys
import sqlite3
import click
from pymlst.lib import sql

desc = "Remove gene to a wgMLST database and sequences specificaly associated"


@click.command()
@click.option('--list', '-l',
              type=click.File('r'), default=None,
              help='File list of genes to removed on the wgMLST database')
@click.argument('database',
                type=click.File('r'), nargs=1)
@click.argument('genes',
                type=str, nargs=-1)
def cli(list, genes, database):
    """Remove gene to a wgMLST database"""

    # list genes to removed

    all_genes = []
    if list is not None:
        for line in list.readlines():
            all_genes.append(line.rstrip("\n"))
    if genes is not None:
        for gene in genes:
            all_genes.append(gene)
    if len(all_genes) == 0:
        raise Exception("No gene to removed found.\n")
    all_genes = set(all_genes)

    try:
        db = sqlite3.connect(database.name)
        cursor = db.cursor()

        # index for old database
        sql.index_database(cursor)

        for gene in all_genes:
            sys.stderr.write(gene + "     ")

            # Search seq ids
            cursor.execute('''SELECT seqid FROM mlst WHERE gene=?''', (gene,))
            seqids = cursor.fetchall()
            if len(seqids) == 0:
                raise Exception("Gene name not found in database\n" + gene)

            # remove sample
            cursor.execute('''DELETE FROM mlst WHERE gene=? ''', (gene,))

            # remove seqs if no other gene have this seq
            for seqid in seqids:
                cursor.execute('''DELETE from sequences as s
                                  where not exists (
                                  select 1 from mlst where seqid=s.id)
                                  and id=?''', (seqid[0],))
            sys.stderr.write("OK\n")

        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
